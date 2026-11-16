"""Reproduce harmonic/phase examples using archived fifth-order gauge histories."""

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from analyze_case import accuracy_metrics
from harmonic_analysis import design_matrix, harmonic_analysis
from reflection.reflection_coeff_komega import load_profiles


ROOT = Path(__file__).resolve().parent
RESOLUTIONS = [
    ("post_processing_m2levs", 0.25),
    ("post_processing_m1levs", 0.125),
    ("post_processing_0levs", 0.0625),
    ("post_processing_1levs", 0.03125),
    ("post_processing_2levs", 0.015625),
]


def load_examples(data, start=10.0, end=18.0):
    options = dict(dt=0.001, output_interval=10, start=start, end=end,
                   x_min=9.99, x_max=10.01, prefix="sampling_fs_gauge",
                   allow_single_point=True)
    x, times, target = load_profiles(data / "ow_vof/post_processing", **options)
    examples = []
    for label, dx in RESOLUTIONS:
        xx, tt, eta = load_profiles(data / "vof" / label, **options)
        if not np.array_equal(tt, times) or not np.allclose(x, xx, rtol=0, atol=1e-8):
            raise ValueError(f"{label}: archived gauge does not match target coordinates/times")
        summary, arrays = harmonic_analysis(times, eta, target)
        examples.append({
            "label": f"H/dx = {0.3585 / dx:.2f}", "dx_m": dx,
            "cells_per_height": 0.3585 / dx,
            "harmonic_analysis": summary, "arrays": arrays, "eta": eta,
            "accuracy": accuracy_metrics(eta, target, 0.3585),
        })
    return times, target, examples


def make_plots(times, target, examples, output):
    import matplotlib.pyplot as plt

    colors = plt.get_cmap("viridis")(np.linspace(0.1, 0.9, len(examples)))
    orders = examples[0]["arrays"]["harmonic_orders"]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), layout="constrained")
    target_rows = examples[0]["harmonic_analysis"]["harmonics"]
    axes[0].semilogy(orders, [row["target_amplitude_rms_m"] for row in target_rows],
                     "ks--", linewidth=2, label="Archived target")
    for color, example in zip(colors, examples):
        rows = example["harmonic_analysis"]["harmonics"]
        axes[0].semilogy(orders, [row["signal_amplitude_rms_m"] for row in rows],
                         "o-", color=color, label=example["label"])
        axes[1].plot(orders, [row["phase_error_circular_mean_deg"] for row in rows],
                     "o-", color=color)
        axes[2].plot(times, example["eta"][:, 0] - target[:, 0], color=color)
    axes[0].set(xlabel="Harmonic order", ylabel="Amplitude [m]", xticks=orders)
    axes[0].legend(fontsize=8)
    axes[1].set(xlabel="Harmonic order", ylabel="Wrapped phase error [deg]; + = lag", xticks=orders)
    axes[2].set(xlabel="Time [s]", ylabel="Computed minus target elevation [m]")
    for ax in axes:
        ax.grid(True, alpha=0.25)
    axes[1].axhline(0, color="gray", linewidth=0.7)
    axes[2].axhline(0, color="gray", linewidth=0.7)
    frequency = examples[0]["harmonic_analysis"]["fundamental_frequency_hz"]
    fig.suptitle(
        f"ARCHIVED fifth-order Stokes: x=10 m, {times[0]:g} <= t < {times[-1] + times[1] - times[0]:g} s\n"
        f"Common target frequency {frequency:.6f} Hz; about {len(times) * (times[1] - times[0]) * frequency:.1f} periods. "
        "Not the new sweep results."
    )
    fig.savefig(output / "archived_harmonics_phase.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    axes[0, 0].plot(times, target[:, 0], "k--", linewidth=2, label="Archived target")
    for index in (0, 2, 4):
        example = examples[index]
        axes[0, 0].plot(times, example["eta"][:, 0], color=colors[index],
                        label=example["label"])
    axes[0, 0].set(xlim=(times[-1] - 3.2, times[-1]),
                      xlabel="Time [s]", ylabel="Surface elevation [m]")
    axes[0, 0].legend(fontsize=8)
    resolution = [example["cells_per_height"] for example in examples]
    for h in range(5):
        amplitude_errors = [
            100 * example["harmonic_analysis"]["harmonics"][h]["relative_amplitude_error"]
            if example["harmonic_analysis"]["harmonics"][h]["relative_amplitude_error"] is not None else np.nan
            for example in examples
        ]
        phase_errors = [
            example["harmonic_analysis"]["harmonics"][h]["phase_error_circular_mean_deg"]
            for example in examples
        ]
        axes[0, 1].semilogx(resolution, amplitude_errors, "o-", label=f"h={h + 1}")
        axes[1, 0].semilogx(resolution, phase_errors, "o-", label=f"h={h + 1}")
    axes[0, 1].set(xlabel="Cells per wave height, H/dx", ylabel="Harmonic amplitude error [%]")
    axes[1, 0].set(xlabel="Cells per wave height, H/dx", ylabel="Wrapped phase error [deg]; + = lag")
    axes[0, 1].legend(ncol=2, fontsize=8)
    axes[1, 0].legend(ncol=2, fontsize=8)
    fine = examples[-1]
    matrix = design_matrix(times, frequency, 5)
    coefficients = fine["arrays"]["harmonic_signal_coefficient_m"][:, 0]
    oscillatory = matrix[:, 1::2] @ coefficients.real + matrix[:, 2::2] @ coefficients.imag
    reconstructed = oscillatory + np.mean(fine["eta"][:, 0] - oscillatory)
    axes[1, 1].plot(times, (fine["eta"][:, 0] - reconstructed) * 1000, color=colors[-1])
    axes[1, 1].set(xlabel="Time [s]", ylabel="Finest-run harmonic-fit residual [mm]")
    for ax in axes.ravel():
        ax.grid(True, alpha=0.25)
    axes[0, 1].axhline(0, color="gray", linewidth=0.7)
    axes[1, 0].axhline(0, color="gray", linewidth=0.7)
    fig.suptitle("ARCHIVED fifth-order comparison: amplitude, phase and residuals\n"
                 "Short record, historical fixed dt; phase includes reflection/propagation and is not aligned away.")
    fig.savefig(output / "archived_harmonics_convergence.png", dpi=180)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=ROOT / "data/convergence")
    parser.add_argument("--output", type=Path, default=ROOT / "plotting_outputs")
    parser.add_argument("--start", type=float, default=10.0)
    parser.add_argument("--end", type=float, default=18.0)
    args = parser.parse_args(argv)
    times, target, examples = load_examples(args.data, args.start, args.end)
    args.output.mkdir(parents=True, exist_ok=True)
    reports = [
        {key: value for key, value in example.items() if key not in ("arrays", "eta")}
        for example in examples
    ]
    report = {
        "description": "Historical fifth-order gauge runs, not newly generated sweep results",
        "source": str(args.data.resolve()), "gauge_x_m": 10.0,
        "time_window_s": [args.start, args.end],
        "reference": "Stored ow_vof gauge, shared by archived mesh comparisons; not assumed same-grid",
        "cases": reports,
    }
    (args.output / "archived_harmonics.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    with (args.output / "archived_harmonics.csv").open("w", newline="") as stream:
        columns = ["cells_per_height", "order", "signal_amplitude_rms_m", "target_amplitude_rms_m",
                   "relative_amplitude_error", "phase_error_circular_mean_deg", "equivalent_time_lag_s",
                   "phase_status"]
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for example in examples:
            for row in example["harmonic_analysis"]["harmonics"]:
                writer.writerow(dict(cells_per_height=example["cells_per_height"], **row))
    np.savez_compressed(
        args.output / "archived_harmonics.npz",
        times_s=times, target_eta_m=target, cells_per_height=[e["cells_per_height"] for e in examples],
        signal_eta_m=np.stack([e["eta"] for e in examples]),
        **{
            key: np.stack([e["arrays"][key] for e in examples])
            for key in examples[0]["arrays"]
        },
    )
    make_plots(times, target, examples, args.output)
    print(f"Saved archived harmonic diagnostics in {args.output}")
    for example in examples:
        fundamental = example["harmonic_analysis"]["harmonics"][0]
        print(f"{example['label']}: A1 error={100 * fundamental['relative_amplitude_error']:.3f}%, "
              f"phase={fundamental['phase_error_circular_mean_deg']:.3f} deg, "
              f"lag={1000 * fundamental['equivalent_time_lag_s']:.3f} ms")


if __name__ == "__main__":
    main()
