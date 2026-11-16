"""Directional k-omega reflection analysis of ASCII FreeSurfaceSampler profiles."""

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import re
import sys

import numpy as np


def uniform_spacing(values, name):
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or values.size < 8 or not np.isfinite(values).all():
        raise ValueError(f"{name} must contain at least eight finite coordinates")
    spacing = np.diff(values)
    if spacing[0] <= 0 or not np.allclose(spacing, spacing[0], rtol=1e-6, atol=1e-10):
        raise ValueError(f"{name} must be increasing and uniformly spaced")
    return float(spacing[0])


def load_profiles(directory, *, dt, output_interval, start, end, x_min, x_max,
                  prefix="sampling_fs_domain", allow_single_point=False):
    """Load eta[t, x]; filenames carry step numbers, not physical timestamps."""
    if not np.isfinite([dt, start, end, x_min, x_max]).all():
        raise ValueError("Sampling parameters must be finite")
    if dt <= 0 or output_interval <= 0 or end <= start or x_max <= x_min:
        raise ValueError("Require positive sampling intervals and increasing time/space bounds")
    snapshots = []
    for path in Path(directory).glob(f"{prefix}*.txt"):
        match = re.fullmatch(re.escape(prefix) + r"(\d+)\.txt", path.name)
        if match:
            step = int(match[1])
            time = step * dt
            if start - 1e-10 <= time < end - 1e-10:
                snapshots.append((step, path))
    snapshots.sort()
    if len(snapshots) < 8:
        raise ValueError(f"{directory}: fewer than eight profiles in [{start}, {end}) s")
    steps = np.array([step for step, _ in snapshots])
    times = steps * dt
    sample_dt = dt * output_interval
    if not np.all(np.diff(steps) == output_interval):
        raise ValueError("Missing/duplicate profiles or output interval does not match filenames")
    if times[0] >= start + sample_dt - 1e-10 or times[-1] < end - sample_dt - 1e-10:
        raise ValueError("Profiles do not cover the requested time window")

    reference = None
    elevations = []
    for _, path in snapshots:
        with path.open() as stream:
            count = int(stream.readline().strip())
        data = np.loadtxt(path, skiprows=5, usecols=(0, 1, 2), ndmin=2)
        if data.shape[0] != count or not np.isfinite(data).all():
            raise ValueError(f"{path}: invalid point count or nonfinite surface coordinates")
        data = data[np.argsort(data[:, 0])]
        if not np.allclose(data[:, 1], data[0, 1], rtol=0, atol=1e-10):
            raise ValueError(f"{path}: expected a single streamwise line, not multiple y values")
        single_point = data.shape[0] == 1
        if single_point and not allow_single_point:
            raise ValueError(f"{path}: single-point gauge requires allow_single_point=True")
        if not single_point:
            dx = uniform_spacing(data[:, 0], str(path))
        if reference is None:
            reference = data[:, :2].copy()
            if not single_point and (data[0, 0] > x_min + 1e-8 or data[-1, 0] < x_max - dx - 1e-8):
                raise ValueError(f"{path}: profiles do not cover the requested spatial window")
            selected = (data[:, 0] >= x_min - 1e-8) & (data[:, 0] < x_max - 1e-8)
            x = data[selected, 0]
            if not len(x):
                raise ValueError(f"{path}: no gauge/profile points in the requested spatial window")
            if not single_point:
                uniform_spacing(x, "selected x")
        elif data[:, :2].shape != reference.shape or not np.allclose(
            data[:, :2], reference, rtol=0, atol=1e-8
        ):
            raise ValueError(f"{path}: surface sampling coordinates changed between snapshots")
        elevations.append(data[selected, 2])
    return x, times, np.stack(elevations)


@dataclass
class Decomposition:
    k: np.ndarray
    omega: np.ndarray
    variance: np.ndarray
    incident_variance: float
    reflected_variance: float

    @property
    def reflection_coefficient(self):
        return float(np.sqrt(self.reflected_variance / self.incident_variance))


def decompose(x, times, eta, frequency_band=(0.4, 0.85), incident_direction=1):
    """Integrate directional elevation variance in a frequency band (Hz).

    Transform convention: exp(-i*k*x + i*omega*t), so k > 0 at omega > 0
    travels toward +x. Zero and spatial Nyquist modes have no assigned direction.
    """
    dx = uniform_spacing(x, "x")
    dt = uniform_spacing(times, "time")
    eta = np.asarray(eta, dtype=float)
    nt, nx = len(times), len(x)
    if eta.shape != (nt, nx) or not np.isfinite(eta).all():
        raise ValueError("eta must be finite with shape (number of times, number of x points)")
    low, high = frequency_band
    if not np.isfinite([low, high]).all() or not 0 < low < high < 0.5 / dt:
        raise ValueError("Frequency band must lie strictly between zero and temporal Nyquist")
    if incident_direction not in (-1, 1):
        raise ValueError("incident_direction must be +1 or -1")
    # Periodic Hann windows; mean-square normalization preserves integrated variance.
    wt = np.hanning(nt + 1)[:-1]
    wx = np.hanning(nx + 1)[:-1]
    centered = eta - eta.mean(axis=0, keepdims=True)
    transformed = np.fft.fft(
        np.fft.ifft(centered * wt[:, None] * wx[None, :], axis=0), axis=1
    ) * nt
    k = 2 * np.pi * np.fft.fftfreq(nx, dx)
    frequency = np.fft.fftfreq(nt, dt)
    positive = frequency > 0
    frequency = frequency[positive]
    variance = 2 * np.abs(transformed[positive]) ** 2 / (
        (nt * nx) ** 2 * np.mean(wt**2) * np.mean(wx**2)
    )
    band = (frequency >= low) & (frequency <= high)
    if not band.any():
        raise ValueError("No FFT frequencies in the selected band; use a longer time window")
    directional = (k != 0) & (np.abs(k) < np.pi / dx * (1 - 1e-10))
    incident = directional & (incident_direction * k > 0)
    reflected = directional & (incident_direction * k < 0)
    ei = float(variance[np.ix_(band, incident)].sum())
    er = float(variance[np.ix_(band, reflected)].sum())
    if ei <= np.finfo(float).tiny:
        raise ValueError("No incident energy in the selected frequency band")
    order = np.argsort(k)
    return Decomposition(k[order], 2 * np.pi * frequency, variance[:, order], ei, er)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="directory containing sampling_fs_domain*.txt")
    parser.add_argument("--dt", type=float, default=0.001, help="simulation fixed dt in seconds")
    parser.add_argument("--output-interval", type=int, default=10, help="sampler interval in steps")
    parser.add_argument("--start", type=float, default=40.0)
    parser.add_argument("--end", type=float, default=72.0, help="exclusive time bound in seconds")
    parser.add_argument("--x-min", type=float, default=4.0)
    parser.add_argument("--x-max", type=float, default=12.0, help="exclusive spatial bound in meters")
    parser.add_argument("--frequency-band", type=float, nargs=2, default=(0.4, 0.85), metavar=("LOW", "HIGH"))
    parser.add_argument("--incident-direction", choices=("+x", "-x"), default="+x")
    parser.add_argument("--output", type=Path, required=True, help="output prefix for .json and .npz")
    args = parser.parse_args(argv)
    x, times, eta = load_profiles(
        args.directory, dt=args.dt, output_interval=args.output_interval,
        start=args.start, end=args.end, x_min=args.x_min, x_max=args.x_max,
    )
    result = decompose(x, times, eta, args.frequency_band, 1 if args.incident_direction == "+x" else -1)
    duration = len(times) * uniform_spacing(times, "time")
    length = len(x) * uniform_spacing(x, "x")
    if duration < 16:
        print("WARNING: fewer than ten nominal wave periods; reflection may be transient "
              "or dominated by spectral leakage.", file=sys.stderr)
    if length < 8 - 1e-8:
        print("WARNING: fewer than two nominal wavelengths; directional leakage may be large.",
              file=sys.stderr)
    summary = {
        "source": str(args.directory.resolve()),
        "time_first_s": float(times[0]), "time_last_s": float(times[-1]),
        "x_first_m": float(x[0]), "x_last_m": float(x[-1]),
        "nt": len(times), "nx": len(x),
        "delta_f_hz": 1 / duration, "delta_k_rad_per_m": 2 * np.pi / length,
        "frequency_band_hz": list(args.frequency_band),
        "incident_direction": args.incident_direction,
        "incident_variance_m2": result.incident_variance,
        "reflected_variance_m2": result.reflected_variance,
        "incident_equivalent_amplitude_m": float(np.sqrt(2 * result.incident_variance)),
        "reflected_equivalent_amplitude_m": float(np.sqrt(2 * result.reflected_variance)),
        "reflection_coefficient": result.reflection_coefficient,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.with_suffix(".json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    np.savez_compressed(
        args.output.with_suffix(".npz"), k_rad_per_m=result.k,
        omega_rad_per_s=result.omega, variance_m2=result.variance,
    )
    print(json.dumps(summary, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
