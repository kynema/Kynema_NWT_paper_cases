"""Largest-wave and crest comparisons between HOS-NWT and Kynema-SGF only."""

import argparse
import csv
from dataclasses import dataclass
import json
from pathlib import Path
import re
import warnings

import numpy as np


ROOT = Path(__file__).resolve().parent
COLORS = {"hos": "black", "sgf_coarse": "#2878b5", "sgf_fine": "#e17c05"}


@dataclass
class Record:
    key: str
    label: str
    time: np.ndarray
    eta: np.ndarray
    metadata: dict


def validate_record(record):
    t, z = record.time, record.eta
    if (t.ndim != 1 or t.shape != z.shape or len(t) < 8
            or not np.isfinite(t).all() or not np.isfinite(z).all()):
        raise ValueError(f"{record.label}: require matching finite time/elevation vectors")
    if np.any(np.diff(t) < 0) or t[-1] <= t[0]:
        raise ValueError(f"{record.label}: timestamps must be nondecreasing and span a positive interval")


def unique_samples(record):
    """Average tied printed timestamps for interpolation, never crest extraction."""
    starts = np.r_[0, np.flatnonzero(np.diff(record.time) > 0) + 1]
    counts = np.diff(np.r_[starts, len(record.time)])
    return record.time[starts], np.add.reduceat(record.eta, starts) / counts


def load_hos(path, probe=2):
    if probe < 1:
        raise ValueError("HOS probe indices start at one")
    header = []
    with Path(path).open() as stream:
        for line_number, line in enumerate(stream):
            stripped = line.strip()
            if not stripped or stripped.startswith(("#", "TITLE", "VARIABLES")):
                header.append(line)
                continue
            first_time = stripped.split()[0]
            break
        else:
            raise ValueError(f"{path}: missing HOS numeric records")
    text = "".join(header)

    def parameter(name):
        match = re.search(rf"\b{re.escape(name)}::\s*([0-9.Ee+-]+)", text)
        if match is None:
            raise ValueError(f"{path}: missing {name} in HOS header")
        value = float(match[1])
        if not np.isfinite(value):
            raise ValueError(f"{path}: nonfinite HOS {name}")
        return value

    frequency, peak_period = parameter("f_out"), parameter("Tp")
    if frequency <= 0 or peak_period <= 0 or parameter("idim") != 1:
        raise ValueError("Require positive HOS output frequency/peak period and dimensional records")
    data = np.loadtxt(path, skiprows=line_number, usecols=(0, probe), ndmin=2)
    printed, eta = data.T
    precision = re.fullmatch(r"[+-]?\d+\.(\d+)[Ee][+-]?\d+", first_time)
    if precision is None:
        raise ValueError("Expected scientific-notation HOS timestamps with explicit printed precision")
    decimals = len(precision[1])
    time = printed[0] + np.arange(len(printed)) / frequency
    tolerance = 0.50001 * 10.0**(np.floor(np.log10(np.maximum(abs(printed), 1e-300))) - decimals)
    if not np.isfinite(printed).all() or np.any(abs(time - printed) > tolerance + 1e-10):
        raise ValueError("HOS timestamps disagree with the header output frequency beyond rounding precision")
    record = Record("hos", f"HOS-NWT precursor (p{probe})", time, eta, {
        "source": str(Path(path).resolve()), "probe": probe, "output_frequency_hz": frequency,
        "peak_period_s": peak_period, "printed_time_decimal_digits": decimals,
        "printed_time_ties": int(np.sum(np.diff(printed) == 0)),
        "clock_policy": "t[0] + sample_index/f_out, checked against printed rounding intervals",
    })
    validate_record(record)
    return record


def load_sgf(path, key, label):
    data = np.loadtxt(path, comments="#", usecols=(0, 1, 2, 3), ndmin=2)
    if not np.isfinite(data).all() or not np.allclose(data[:, 1:3], [0, 0], rtol=0, atol=1e-10):
        raise ValueError(f"{path}: expected a finite, fixed x=y=0 gauge")
    record = Record(key, label, data[:, 0], data[:, 3], {
        "source": str(Path(path).resolve()), "gauge_xy_m": [0.0, 0.0],
        "printed_time_ties": int(np.sum(np.diff(data[:, 0]) == 0)),
        "clock_policy": "native printed times; tied elevations retained for crest extraction",
    })
    validate_record(record)
    return record


def clip_record(record, start, end):
    t, z = unique_samples(record)
    if not np.isfinite([start, end]).all() or not t[0] <= start < end <= t[-1]:
        raise ValueError(f"{record.label}: requested interval is outside available data")
    interior = (record.time > start) & (record.time < end)
    clipped = Record(
        record.key, record.label,
        np.r_[start, record.time[interior], end],
        np.r_[np.interp(start, t, z), record.eta[interior], np.interp(end, t, z)],
        record.metadata,
    )
    validate_record(clipped)
    return clipped


def water_level_moments(record):
    t, z = unique_samples(record)
    dt = np.diff(t)
    mean = float(np.sum(dt * (z[:-1] + z[1:]) / 2) / (t[-1] - t[0]))
    a, b = z[:-1] - mean, z[1:] - mean
    variance = float(np.sum(dt * (a*a + a*b + b*b) / 3) / (t[-1] - t[0]))
    return mean, np.sqrt(variance)


def upcrossing_waves(record, mean):
    """One native-sample crest/trough per complete mean-level upcrossing wave."""
    centered = record.eta - mean
    indices = np.flatnonzero((centered[:-1] <= 0) & (centered[1:] > 0))
    crossing_times = record.time[indices] - centered[indices] * (
        record.time[indices + 1] - record.time[indices]
    ) / (centered[indices + 1] - centered[indices])
    if np.any(np.diff(crossing_times) <= 0):
        raise ValueError(f"{record.label}: multiple waves within one printed timestamp; clock resolution is inadequate")
    waves = []
    for i, (left, right) in enumerate(zip(indices[:-1], indices[1:])):
        segment = centered[left + 1:right + 1]
        crest = left + 1 + int(np.argmax(segment))
        trough = left + 1 + int(np.argmin(segment))
        crest_height, trough_height = float(centered[crest]), min(0.0, float(centered[trough]))
        waves.append({
            "wave_index": i, "start_s": float(crossing_times[i]),
            "end_s": float(crossing_times[i + 1]),
            "period_s": float(crossing_times[i + 1] - crossing_times[i]),
            "crest_time_s": float(record.time[crest]), "crest_height_m": crest_height,
            "crest_elevation_m": float(record.eta[crest]),
            "trough_time_s": float(record.time[trough]), "trough_height_m": trough_height,
            "wave_height_m": crest_height - trough_height,
        })
    if not waves:
        raise ValueError(f"{record.label}: no complete mean-level upcrossing waves")
    return waves


def select_events(waves, start, end, half_window, count, rank_by="crest"):
    if not np.isfinite(half_window) or half_window <= 0 or count < 1:
        raise ValueError("Require positive event count and half-window")
    if rank_by not in ("crest", "height"):
        raise ValueError("Rank events by crest or height")
    field = "crest_height_m" if rank_by == "crest" else "wave_height_m"
    candidates = sorted(waves, key=lambda wave: wave[field], reverse=True)
    selected = []
    for wave in candidates:
        center = wave["crest_time_s"]
        if center - half_window < start or center + half_window > end:
            continue
        if all(abs(center - other["crest_time_s"]) >= 2 * half_window for other in selected):
            selected.append(wave)
        if len(selected) == count:
            break
    if not selected:
        raise ValueError("No complete HOS events with full surrounding context")
    if len(selected) != count:
        warnings.warn(f"Only {len(selected)} nonoverlapping HOS events available, requested {count}", RuntimeWarning)
    return selected


def agreement(signal, reference):
    if signal.shape != reference.shape or not len(signal):
        raise ValueError("Temporal comparison requires matching nonempty histories")
    a, b = signal - signal.mean(), reference - reference.mean()
    denominator = np.sqrt(np.mean(a*a) * np.mean(b*b))
    if denominator <= 0:
        raise ValueError("Temporal correlation is undefined for a constant history")
    difference = signal - reference
    return {
        "bias_m": float(np.mean(difference)), "rmse_m": float(np.sqrt(np.mean(difference**2))),
        "mae_m": float(np.mean(abs(difference))),
        "correlation": float(np.mean(a*b) / denominator),
    }


def compare(records, *, start, end, peak_period, event_count=6, half_window_periods=3,
            sample_dt=0.01, rank_by="crest"):
    if (not np.isfinite([start, end, peak_period, sample_dt, half_window_periods]).all()
            or peak_period <= 0 or sample_dt <= 0 or end <= start):
        raise ValueError("Require finite increasing time bounds and positive period/sample interval")
    if [record.key for record in records] != ["hos", "sgf_coarse", "sgf_fine"]:
        raise ValueError("Expected HOS-NWT, coarse SGF and fine SGF in that order")
    records = [clip_record(record, start, end) for record in records]
    time = start + np.arange(int(np.floor((end - start) / sample_dt + 1e-9)) + 1) * sample_dt
    if len(time) < 8:
        raise ValueError("Comparison time grid has fewer than eight samples")
    time = time[time <= end]
    arrays = {"comparison_time_s": time}
    statistics, all_waves = {}, {}
    for record in records:
        mean, sigma = water_level_moments(record)
        waves = upcrossing_waves(record, mean)
        crests = np.array([wave["crest_height_m"] for wave in waves])
        all_waves[record.key] = waves
        unique_t, unique_eta = unique_samples(record)
        arrays[f"{record.key}_eta_m"] = np.interp(time, unique_t, unique_eta)
        arrays[f"{record.key}_crest_height_m"] = crests
        arrays[f"{record.key}_crest_time_s"] = np.array([wave["crest_time_s"] for wave in waves])
        statistics[record.key] = {
            "label": record.label, "complete_wave_count": len(waves),
            "mean_water_level_m": mean, "surface_std_m": float(sigma), "hm0_m": float(4 * sigma),
            "mean_crest_height_m": float(crests.mean()),
            **{f"crest_p{label}_m": float(np.quantile(crests, p))
               for label, p in (("50", 0.5), ("90", 0.9), ("95", 0.95), ("99", 0.99), ("99p9", 0.999))},
            "maximum_crest_height_m": float(crests.max()),
            "maximum_wave_height_m": float(max(wave["wave_height_m"] for wave in waves)),
            "native_maximum_elevation_m": float(record.eta.max()),
            "native_maximum_time_s": float(record.time[record.eta.argmax()]),
            "tied_samples_in_window": int(np.sum(np.diff(record.time) == 0)),
            "minimum_wave_period_s": float(min(wave["period_s"] for wave in waves)),
            "waves_shorter_than_half_peak_period": sum(wave["period_s"] < peak_period / 2 for wave in waves),
        }
        if record.key != "hos":
            statistics[record.key]["temporal_agreement"] = agreement(
                arrays[f"{record.key}_eta_m"], arrays["hos_eta_m"]
            )
    half_window = half_window_periods * peak_period
    chosen = select_events(all_waves["hos"], start, end, half_window, event_count, rank_by)
    events = []
    for rank, wave in enumerate(chosen, 1):
        center = wave["crest_time_s"]
        context = (time >= center - half_window) & (time <= center + half_window)
        event = {"rank": rank, "hos_wave": wave, "series": {}}
        for record in records:
            mask = (record.time >= wave["start_s"]) & (record.time <= wave["end_s"])
            indices = np.flatnonzero(mask)
            if len(indices) < 3:
                raise ValueError(f"{record.label}: insufficient native samples in selected HOS wave")
            local = int(np.argmax(record.eta[indices]))
            peak = indices[local]
            height = float(record.eta[peak] - statistics[record.key]["mean_water_level_m"])
            event["series"][record.key] = {
                "crest_time_s": float(record.time[peak]),
                "crest_elevation_m": float(record.eta[peak]),
                "crest_height_m": height,
                "crest_height_error_m": height - wave["crest_height_m"],
                "crest_height_relative_error": height / wave["crest_height_m"] - 1,
                "crest_time_lag_s": float(record.time[peak] - center),
                "peak_at_matching_interval_edge": local in (0, len(indices) - 1),
                "context_agreement": agreement(
                    arrays[f"{record.key}_eta_m"][context], arrays["hos_eta_m"][context]
                ),
            }
        events.append(event)
    report = {
        "analysis_window_s": [start, end], "peak_period_s": peak_period,
        "comparison_sample_dt_s": sample_dt, "event_half_window_s": half_window,
        "event_rank_by": rank_by,
        "event_selection": f"largest HOS {'mean-referenced crests' if rank_by == 'crest' else 'upcrossing wave heights'}, nonoverlapping full-context windows",
        "crest_definition": "native-sample maximum above each record's time-weighted mean in each complete upcrossing wave",
        "clock_alignment": "original physical clock, no fitted phase/time shifts",
        "timestamp_ties": "native elevations preserved for crests; tied times averaged only for interpolation/moments",
        "gauge_mapping": "HOS p2 follows the existing Julia analysis; expected x=16 m maps to SGF x=0 via offset -16 m. The archived probe_2.inp is not present; verify that position list independently.",
        "sources": {record.key: record.metadata for record in records},
        "statistics": statistics, "events": events,
    }
    return report, arrays, records, all_waves


def make_plots(report, arrays, records, output):
    import matplotlib.pyplot as plt

    events = report["events"]
    rows = (len(events) + 1) // 2
    fig, axes = plt.subplots(rows, 2, figsize=(13, 3.1 * rows + 1), squeeze=False, sharex=True)
    half_window = report["event_half_window_s"]
    for ax, event in zip(axes.ravel(), events):
        wave = event["hos_wave"]
        center = wave["crest_time_s"]
        for record in records:
            mask = abs(record.time - center) <= half_window
            ax.plot(record.time[mask] - center, record.eta[mask], color=COLORS[record.key],
                    linewidth=1.3, linestyle="--" if record.key == "sgf_coarse" else "-", label=record.label)
            peak = event["series"][record.key]
            ax.plot(peak["crest_time_s"] - center, peak["crest_elevation_m"], "o",
                    color=COLORS[record.key], markersize=4)
        ax.axvspan(wave["start_s"] - center, wave["end_s"] - center, color="gray", alpha=0.1)
        ax.axvline(0, color="gray", linewidth=0.6)
        ax.set(title=f"#{event['rank']}: HOS crest at {center:.3f} s, C={wave['crest_height_m']:.3f} m",
               ylabel="Elevation above still water [m]", xlim=(-half_window, half_window))
        ax.grid(True, alpha=0.2)
    for ax in axes.ravel()[len(events):]:
        ax.set_visible(False)
    for ax in axes[-1]:
        ax.set_xlabel("Time relative to HOS crest [s]; no phase alignment")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.955), ncol=3, fontsize=9)
    metric = "crest height" if report["event_rank_by"] == "crest" else "wave height"
    fig.suptitle(f"Irregular tank: separated HOS-NWT events ranked by {metric}", y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(output / "largest_wave_events.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    reference = arrays["hos_crest_height_m"]
    probabilities = np.linspace(0.5, 0.999, 100)
    reference_q = np.quantile(reference, probabilities)
    groups = np.arange(4)
    for i, record in enumerate(records):
        crests = np.sort(arrays[f"{record.key}_crest_height_m"])
        exceedance = (len(crests) - np.arange(len(crests))) / (len(crests) + 1)
        axes[0, 0].semilogy(crests, exceedance, color=COLORS[record.key], label=record.label)
        if record.key != "hos":
            axes[0, 1].plot(reference_q, np.quantile(crests, probabilities), color=COLORS[record.key],
                            label=record.label)
        values = np.r_[np.quantile(crests, [0.5, 0.9, 0.99]), crests.max()]
        axes[1, 0].bar(groups + (i - 1) * 0.24, values, width=0.23, color=COLORS[record.key],
                       label=record.label)
    axes[0, 1].plot([0, max(reference_q) * 1.05], [0, max(reference_q) * 1.05], "k--", linewidth=1)
    axes[0, 0].set(xlabel="Crest height above each record's mean [m]", ylabel="Crest exceedance (Weibull plotting position)")
    axes[0, 1].set(xlabel="HOS-NWT crest quantile [m]", ylabel="Kynema-SGF crest quantile [m]",
                  title="Distribution Q-Q: 50th to 99.9th percentiles")
    axes[1, 0].set(xticks=groups, xticklabels=["Median", "P90", "P99", "Maximum"],
                  ylabel="Crest height above each record's mean [m]")
    stat = report["statistics"]
    table_rows = [
        ["Complete waves"] + [str(stat[r.key]["complete_wave_count"]) for r in records],
        ["Mean level [mm]"] + [f"{1000 * stat[r.key]['mean_water_level_m']:.2f}" for r in records],
        ["Hm0 = 4 sigma [m]"] + [f"{stat[r.key]['hm0_m']:.4f}" for r in records],
        ["P99 crest [m]"] + [f"{stat[r.key]['crest_p99_m']:.4f}" for r in records],
        ["Max crest [m]"] + [f"{stat[r.key]['maximum_crest_height_m']:.4f}" for r in records],
    ]
    axes[1, 1].axis("off")
    table = axes[1, 1].table(cellText=table_rows, colLabels=["Metric", "HOS", "SGF coarse", "SGF fine"],
                              loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.7)
    axes[1, 1].set_title("Native crests; shared record; no experiment")
    for ax in (axes[0, 0], axes[0, 1], axes[1, 0]):
        ax.grid(True, alpha=0.2)
    axes[0, 0].legend(fontsize=8)
    axes[0, 1].legend(fontsize=8)
    start, end = report["analysis_window_s"]
    fig.suptitle(f"Irregular-tank crest statistics: {start:g}--{end:g} s\n"
                 "One crest per complete mean-upcrossing wave; extreme tail is sample-limited.")
    fig.savefig(output / "crest_height_statistics.png", dpi=180)
    plt.close(fig)


def save_results(report, arrays, records, waves, output, plots=True):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    (output / "temporal_comparison.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    np.savez_compressed(output / "temporal_comparison.npz", **arrays)
    statistics = []
    event_rows = []
    crest_rows = []
    for record in records:
        stats = report["statistics"][record.key]
        statistics.append(dict(source=record.key, **{key: value for key, value in stats.items()
                                                    if key != "temporal_agreement"},
                               **stats.get("temporal_agreement", {})))
        crest_rows.extend(dict(source=record.key, **wave) for wave in waves[record.key])
        for event in report["events"]:
            peak = event["series"][record.key]
            event_rows.append(dict(
                event_rank=event["rank"], source=record.key,
                hos_crest_time_s=event["hos_wave"]["crest_time_s"],
                **{key: value for key, value in peak.items() if key != "context_agreement"},
                **peak["context_agreement"],
            ))
    for name, rows in (("crest_statistics.csv", statistics), ("largest_wave_events.csv", event_rows),
                       ("wave_crests.csv", crest_rows)):
        columns = list(dict.fromkeys(key for row in rows for key in row))
        with (output / name).open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    if plots:
        make_plots(report, arrays, records, output)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=ROOT / "data")
    parser.add_argument("--output", type=Path, default=ROOT / "plotting_outputs")
    parser.add_argument("--start", type=float, default=30.0, help="exclude startup before this model-scale time")
    parser.add_argument("--end", type=float, help="default: end of the common available record")
    parser.add_argument("--events", type=int, default=6)
    parser.add_argument("--rank-by", choices=("crest", "height"), default="crest")
    parser.add_argument("--half-window-periods", type=float, default=3)
    parser.add_argument("--sample-dt", type=float, default=0.01, help="interpolation grid only, not crest extraction")
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args(argv)
    hos = load_hos(args.data / "HOS_NWT_probes.dat")
    records = [
        hos,
        load_sgf(args.data / "amr_wind_hos_nwt_sampling_time_series.txt",
                 "sgf_coarse", "Kynema-SGF, dz=0.0188 m"),
        load_sgf(args.data / "amr_wind_hos_nwt_sampling_time_series_finer.txt",
                 "sgf_fine", "Kynema-SGF, dz=0.0094 m"),
    ]
    end = args.end if args.end is not None else min(record.time[-1] for record in records)
    report, arrays, clipped, waves = compare(
        records, start=args.start, end=end, peak_period=hos.metadata["peak_period_s"],
        event_count=args.events, half_window_periods=args.half_window_periods,
        sample_dt=args.sample_dt, rank_by=args.rank_by,
    )
    save_results(report, arrays, clipped, waves, args.output, plots=not args.no_plots)
    print(f"Compared HOS-NWT and Kynema-SGF only over {args.start:g}--{end:g} s; no clock shifts.")
    for key, stats in report["statistics"].items():
        print(f"{key}: waves={stats['complete_wave_count']}, mean={1000 * stats['mean_water_level_m']:.3f} mm, "
              f"P99 crest={stats['crest_p99_m']:.5f} m, max crest={stats['maximum_crest_height_m']:.5f} m, "
              f"tied samples={stats['tied_samples_in_window']}")
    print(f"Saved {len(report['events'])} event comparisons and crest statistics in {args.output}")


if __name__ == "__main__":
    main()
