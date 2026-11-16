"""Analyze both accuracy and k-omega reflection for a generated regular-wave run."""

import argparse
import json
from pathlib import Path

import numpy as np

from harmonic_analysis import harmonic_analysis, plot_harmonics
from reflection.reflection_coeff_komega import decompose, load_profiles


def accuracy_metrics(eta, target, height):
    eta, target = np.asarray(eta), np.asarray(target)
    if eta.shape != target.shape or not np.isfinite(eta).all() or not np.isfinite(target).all():
        raise ValueError("Surface and target must have matching shapes and finite elevations")
    if not np.isfinite(height) or height <= 0 or not eta.size:
        raise ValueError("Require a positive wave height and nonempty surface histories")
    error = eta - target
    l1 = float(np.mean(np.abs(error)))
    l2 = float(np.sqrt(np.mean(error**2)))
    return {
        "bias_m": float(np.mean(error)),
        "l1_m": l1,
        "l2_m": l2,
        "linf_m": float(np.max(np.abs(error))),
        "l1_over_height": l1 / height,
        "l2_over_height": l2 / height,
        "linf_over_height": float(np.max(np.abs(error))) / height,
    }


def analyze(directory, *, start=None, end=None, harmonic_count=5, phase_floor_ratio=1e-3):
    directory = Path(directory)
    metadata = json.loads((directory / "case.json").read_text())
    start = metadata["analysis_time_s"][0] if start is None else start
    end = metadata["analysis_time_s"][1] if end is None else end
    options = dict(
        dt=metadata["dt_s"], output_interval=metadata["output_interval_steps"],
        start=start, end=end,
        x_min=metadata["sample_x_m"][0], x_max=metadata["sample_x_m"][1],
    )
    x, times, eta = load_profiles(directory / "post_processing", **options)
    xt, tt, target = load_profiles(
        directory / "post_processing", prefix="sampling_fs_target", **options
    )
    if not np.array_equal(times, tt) or x.shape != xt.shape or not np.allclose(x, xt, rtol=0, atol=1e-8):
        raise ValueError("Surface and target coordinates/times do not match")
    spectrum = decompose(x, times, eta, metadata["frequency_band_hz"])
    reference_spectrum = decompose(x, times, target, metadata["frequency_band_hz"])
    harmonics, harmonic_arrays = harmonic_analysis(
        times, eta, target, metadata["frequency_band_hz"],
        count=harmonic_count, phase_floor_ratio=phase_floor_ratio,
    )
    summary = {
        "case": metadata["id"],
        "sweeps": metadata["sweeps"],
        "cells_per_height": metadata["cells_per_height"],
        "aspect_ratio": metadata["aspect_ratio"],
        "nominal_advective_cfl": metadata["nominal_advective_cfl"],
        "dt_s": metadata["dt_s"],
        "time_window_s": [start, end],
        "actual_first_last_time_s": [float(times[0]), float(times[-1])],
        "space_window_m": metadata["sample_x_m"],
        "nt": len(times), "nx": len(x),
        "frequency_band_hz": metadata["frequency_band_hz"],
        "reflection_coefficient": spectrum.reflection_coefficient,
        "incident_variance_m2": spectrum.incident_variance,
        "reflected_variance_m2": spectrum.reflected_variance,
        "target_apparent_reflection_coefficient": reference_spectrum.reflection_coefficient,
        "accuracy_reference": "same-grid solver-prescribed ow_vof surface; phase error retained",
        "accuracy": accuracy_metrics(eta, target, metadata["wave_height_m"]),
        "harmonic_analysis": harmonics,
    }
    return summary, {
        "k_rad_per_m": spectrum.k,
        "omega_rad_per_s": spectrum.omega,
        "variance_m2": spectrum.variance,
        "target_variance_m2": reference_spectrum.variance,
        "x_m": x,
        "rmse_by_x_m": np.sqrt(np.mean((eta - target)**2, axis=0)),
        **harmonic_arrays,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="run directory with case.json and post_processing/")
    parser.add_argument("--start", type=float)
    parser.add_argument("--end", type=float)
    parser.add_argument("--harmonics", type=int, default=5, help="number of jointly fitted harmonics")
    parser.add_argument("--phase-floor-ratio", type=float, default=1e-3)
    parser.add_argument("--plot", action="store_true", help="also save a harmonic/phase PNG (requires Matplotlib)")
    parser.add_argument("--output", type=Path, help="output prefix (default: RUN/analysis/results)")
    args = parser.parse_args(argv)
    summary, arrays = analyze(
        args.directory, start=args.start, end=args.end,
        harmonic_count=args.harmonics, phase_floor_ratio=args.phase_floor_ratio,
    )
    output = args.output if args.output is not None else args.directory / "analysis" / "results"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.with_suffix(".json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    np.savez_compressed(output.with_suffix(".npz"), **arrays)
    if args.plot:
        plot_harmonics(summary, arrays, output.with_name(output.name + "_harmonics").with_suffix(".png"))
    print(json.dumps(summary, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
