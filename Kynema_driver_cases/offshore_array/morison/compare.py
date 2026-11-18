"""Time-aligned engineering/CFD comparisons without fitted phase shifts."""

import csv
import json
from pathlib import Path

import numpy as np

from model import CHANNELS, ORIGINS, array_loads, clock, finite


def read_cfd(path, file_format, sign):
    if sign not in (-1, 1):
        raise ValueError("CFD force-on-structure sign must be +1 or -1")
    if file_format == "csv":
        data = np.genfromtxt(path, delimiter=",", names=True, dtype=float)
        if data.dtype.names != ("time",) + CHANNELS:
            raise ValueError(f"{path}: expected CSV header time," + ",".join(CHANNELS))
        data = np.atleast_1d(data)
        time = data["time"]
        loads = np.column_stack([data[name] for name in CHANNELS])
    elif file_format == "nalu10":
        rows, started = [], False
        with Path(path).open() as stream:
            for number, line in enumerate(stream, 1):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                tokens = line.split()
                try:
                    float(tokens[0])
                except ValueError:
                    if started:
                        raise ValueError(f"{path}:{number}: nonnumeric row inside load history") from None
                    continue
                try:
                    row = [float(token) for token in tokens]
                except ValueError:
                    raise ValueError(f"{path}:{number}: invalid numeric load record") from None
                started = True
                if len(row) != 10:
                    raise ValueError(f"{path}:{number}: nalu10 requires exactly 10 columns")
                rows.append(row)
        if not rows:
            raise ValueError(f"{path}: no numeric load records")
        data = np.array(rows)
        time = data[:, 0]
        loads = np.column_stack((data[:, 1:4] + data[:, 4:7], data[:, 7:10]))
    else:
        raise ValueError(f"Unknown CFD format: {file_format}")
    return clock(time), sign * finite("CFD loads", loads, (len(time), 6))


def welch(reference, cfd, dt):
    """One-sided Hann spectra; cross phase is CFD relative to reference."""
    n = min(4096, len(reference) // 2)
    if n < 8:
        raise ValueError("Need at least 16 common samples for spectral comparison")
    window = np.hanning(n)
    starts = range(0, len(reference) - n + 1, n // 2)
    auto_reference, auto_cfd, cross = [], [], []
    for start in starts:
        a = reference[start:start + n]
        b = cfd[start:start + n]
        fa = np.fft.rfft((a - a.mean(axis=0)) * window[:, None], axis=0)
        fb = np.fft.rfft((b - b.mean(axis=0)) * window[:, None], axis=0)
        auto_reference.append(abs(fa)**2)
        auto_cfd.append(abs(fb)**2)
        cross.append(np.conj(fa) * fb)
    scale = np.full(n // 2 + 1, 2 * dt / np.sum(window**2))
    scale[0] /= 2
    if n % 2 == 0:
        scale[-1] /= 2
    pa = np.mean(auto_reference, axis=0) * scale[:, None]
    pb = np.mean(auto_cfd, axis=0) * scale[:, None]
    ab = np.mean(cross, axis=0) * scale[:, None]
    denominator = pa * pb
    coherence = np.full_like(pa, np.nan)
    np.divide(abs(ab)**2, denominator, out=coherence, where=denominator > 0)
    phase = np.where(denominator > 0, np.rad2deg(np.angle(ab)), np.nan)
    return dict(frequency=np.fft.rfftfreq(n, dt), reference_psd=pa, cfd_psd=pb,
                phase_deg=phase, coherence=np.clip(coherence, 0, 1),
                segment_samples=n, segment_count=len(auto_reference))


def compare(model_path, cfd_paths, output, file_format, sign, start=240, dt=0.15, plot=False):
    if not np.isfinite(dt) or dt <= 0 or not np.isfinite(start):
        raise ValueError("Comparison dt must be positive and start finite")
    if set(cfd_paths) != {1, 2, 3, 4}:
        raise ValueError("Provide one CFD history for each jacket 1 through 4")
    with np.load(model_path, allow_pickle=False) as data:
        time = clock(data["time"])
        loads = finite("model loads", data["loads"], (len(time), 4, 6))
        origins = finite("origins", data["origins"], (4, 3))
        metadata = json.loads(str(data["metadata"]))
    if not np.allclose(origins, ORIGINS, rtol=0, atol=1e-10):
        raise ValueError("Model moments are not about the four configured CFD base origins")
    records = [read_cfd(cfd_paths[jacket], file_format, sign) for jacket in range(1, 5)]
    lower = max([start, time[0]] + [record[0][0] for record in records])
    upper = min([time[-1]] + [record[0][-1] for record in records])
    if upper <= lower:
        raise ValueError("No common analysis interval after startup exclusion")
    common = lower + np.arange(int(np.floor((upper - lower) / dt)) + 1) * dt
    if len(common) < 16:
        raise ValueError("Common analysis interval contains fewer than 16 samples")
    reference = np.empty((len(common), 4, 6))
    cfd = np.empty_like(reference)
    for jacket, (tcfd, lcfd) in enumerate(records):
        for channel in range(6):
            reference[:, jacket, channel] = np.interp(common, time, loads[:, jacket, channel])
            cfd[:, jacket, channel] = np.interp(common, tcfd, lcfd[:, channel])
    reference = np.concatenate((reference, array_loads(reference, origins)[:, None, :]), axis=1)
    cfd = np.concatenate((cfd, array_loads(cfd, origins)[:, None, :]), axis=1)
    labels = ["jacket1", "jacket2", "jacket3", "jacket4", "array_at_world_origin"]
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    statistics, events = [], []
    for jacket, label in enumerate(labels):
        for channel, name in enumerate(CHANNELS):
            a, b = reference[:, jacket, channel], cfd[:, jacket, channel]
            ac, bc = a - a.mean(), b - b.mean()
            denominator = np.linalg.norm(ac) * np.linalg.norm(bc)
            row = dict(body=label, channel=name, units="N" if channel < 3 else "N m",
                       model_mean=float(a.mean()), cfd_mean=float(b.mean()),
                       model_std=float(a.std()), cfd_std=float(b.std()),
                       model_min=float(a.min()), model_max=float(a.max()),
                       cfd_min=float(b.min()), cfd_max=float(b.max()),
                       bias=float((b - a).mean()), rmse=float(np.sqrt(np.mean((b - a)**2))),
                       demeaned_rmse=float(np.sqrt(np.mean((bc - ac)**2))),
                       correlation=float(np.sum(ac * bc) / denominator) if denominator > 0 else None)
            statistics.append(row)
        magnitude = np.linalg.norm(cfd[:, jacket, :2], axis=1)
        candidates = np.flatnonzero((magnitude[1:-1] >= magnitude[:-2]) &
                                   (magnitude[1:-1] > magnitude[2:])) + 1
        selected = []
        for index in candidates[np.argsort(magnitude[candidates])[::-1]]:
            if all(abs(common[index] - other) >= 40 for other in selected):
                selected.append(common[index])
                events.append(dict(body=label, time=float(common[index]),
                                   cfd_horizontal_force=float(magnitude[index]),
                                   model_horizontal_force=float(np.linalg.norm(reference[index, jacket, :2])),
                                   window_start=float(max(lower, common[index] - 20)),
                                   window_end=float(min(upper, common[index] + 20))))
                if len(selected) == 3:
                    break
    for name, rows in [("statistics.csv", statistics), ("largest_load_events.csv", events)]:
        if rows:
            with (output / name).open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
    spectra = welch(reference.reshape(len(common), -1), cfd.reshape(len(common), -1), dt)
    np.savez_compressed(output / "comparison.npz", time=common, reference=reference, cfd=cfd,
                        labels=labels, **spectra)
    report = dict(model_metadata=metadata, interval=[float(common[0]), float(common[-1])],
                  dt=dt, cfd_format=file_format, cfd_sign=sign,
                  cfd_paths={str(key): str(value) for key, value in cfd_paths.items()},
                  moment_origins=origins.tolist(), array_moment_origin=[0, 0, 0],
                  spectral_segment_samples=spectra["segment_samples"],
                  spectral_segment_count=spectra["segment_count"],
                  phase_convention="angle(conj(model FFT) * CFD FFT); positive = CFD leads",
                  warning=metadata.get("scope_warning",
                          "Normal-strip loads only; no hydrostatics or axial/end loads. "
                          "Demeaning does not correct time-varying buoyancy.") + " "
                          "Verify rebuilt CFD geometry, pressure convention, sign and HOS registration. "
                          "Uniform linear interpolation, no anti-alias filter or time-shift fitting.",
                  statistics=statistics, largest_load_events=events)
    (output / "comparison.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    if plot:
        plot_comparison(output, common, reference, cfd, spectra, labels, events, metadata)
    return report


def plot_comparison(output, time, reference, cfd, spectra, labels, events, metadata):
    import matplotlib.pyplot as plt

    reference_label = metadata.get("reference_label", "Morison")
    title = "SYNTHETIC SMOKE TEST" if metadata["source"] == "synthetic" else f"{reference_label} / CFD"
    fig, axes = plt.subplots(4, 2, figsize=(12, 10), constrained_layout=True)
    for jacket in range(4):
        for col, channel in enumerate((0, 4)):
            ax = axes[jacket, col]
            ax.plot(time, reference[:, jacket, channel], label=reference_label)
            ax.plot(time, cfd[:, jacket, channel], label="CFD input", alpha=0.7)
            ax.set_title(f"{labels[jacket]} {CHANNELS[channel]}")
            ax.set_xlabel("Time [s]")
            ax.set_ylabel("N" if channel < 3 else "N m")
    axes[0, 0].legend()
    fig.suptitle(title)
    fig.savefig(output / "time_histories.png", dpi=150)
    plt.close(fig)
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
    for jacket, ax in enumerate(axes.flat):
        index = jacket * 6
        ax.semilogy(spectra["frequency"], spectra["reference_psd"][:, index], label=f"{reference_label} Fx")
        ax.semilogy(spectra["frequency"], spectra["cfd_psd"][:, index], label="CFD input Fx")
        ax.set(xlim=(0, 0.4), xlabel="Frequency [Hz]", ylabel="PSD [N^2/Hz]", title=labels[jacket])
    axes[0, 0].legend()
    fig.suptitle(title)
    fig.savefig(output / "force_spectra.png", dpi=150)
    plt.close(fig)
    selected = [event for event in events if event["body"] == "jacket1"]
    if selected:
        fig, axes = plt.subplots(len(selected), 1, figsize=(10, 3 * len(selected)),
                                 constrained_layout=True, squeeze=False)
        for ax, event in zip(axes[:, 0], selected):
            mask = (time >= event["window_start"]) & (time <= event["window_end"])
            ax.plot(time[mask], reference[mask, 0, 0], label=f"{reference_label} Fx")
            ax.plot(time[mask], cfd[mask, 0, 0], label="CFD input Fx")
            ax.set(xlabel="Time [s]", ylabel="Force [N]", title=f"Jacket 1 event at {event['time']:.2f} s")
        axes[0, 0].legend()
        fig.suptitle(title)
        fig.savefig(output / "largest_load_windows.png", dpi=150)
        plt.close(fig)
