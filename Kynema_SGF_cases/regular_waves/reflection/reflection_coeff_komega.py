"""Directional k-omega reflection analysis of FreeSurfaceSampler profiles (ASCII or native)."""

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


def read_ascii_profile(path):
    """Return (count, x/y/z array, None) from a five-line-header ASCII sampler file.

    Legacy ASCII profiles carry no physical time; callers infer it from the step.
    """
    with path.open() as stream:
        count = int(stream.readline().strip())
    return count, np.loadtxt(path, skiprows=5, usecols=(0, 1, 2), ndmin=2), None


def read_native_time(output_dir):
    """Physical time of a native sampler output from its plotfile-style Header.

    The Header lists a version line, the variable count and names, the
    dimension, then the time at full precision (sampling_info.yaml rounds it
    to six significant digits). Returns None if no Header is present.
    """
    header = output_dir / "Header"
    if not header.is_file():
        return None
    lines = header.read_text().splitlines()
    nvars = int(lines[1])
    return float(lines[2 + nvars + 1])


def read_native_profile(path):
    """Return (count, x/y/z array) from an AMReX native particle sampler directory.

    `path` is the `<label><step>/particles` directory; the output time is read
    from the parent directory's Header when present. The binary layout follows
    kynema-sgf/tools/amrex_particle.py: per grid, (num_ints + 2) int32 values per
    particle, followed by (ndim + num_reals) float64 values per particle whose
    first ndim entries are the particle (sample) position.
    """
    with (path / "Header").open() as fh:
        fh.readline()
        ndim = int(fh.readline())
        num_reals = int(fh.readline())
        for _ in range(num_reals):
            fh.readline()
        num_ints = int(fh.readline())
        for _ in range(num_ints):
            fh.readline()
        fh.readline()  # checkpoint flag
        count = int(fh.readline())
        fh.readline()  # next id
        finest = int(fh.readline())
        ngrids = [int(fh.readline()) for _ in range(finest + 1)]
        grids = [[tuple(int(v) for v in fh.readline().split()) for _ in range(n)] for n in ngrids]
    nints, nreals = num_ints + 2, ndim + num_reals
    rows = []
    for lev, ginfo in enumerate(grids):
        for idx, npts, offset in ginfo:
            if npts < 1:
                continue
            with (path / f"Level_{lev}" / f"DATA_{idx:05d}").open("rb") as fh:
                fh.seek(offset)
                np.fromfile(fh, dtype=np.int32, count=nints * npts)
                reals = np.fromfile(fh, dtype="<f8", count=nreals * npts)
            if reals.size != nreals * npts:
                raise ValueError(f"{path}: truncated particle data in level {lev} grid {idx}")
            rows.append(reals.reshape(npts, nreals)[:, :3])
    data = np.concatenate(rows) if rows else np.empty((0, 3))
    return count, data, read_native_time(path.parent)


def find_profiles(directory, prefix):
    """Map step -> (reader, path) for ASCII files or native particle directories."""
    found = {}
    for path in Path(directory).glob(f"{prefix}*"):
        match = re.fullmatch(re.escape(prefix) + r"(\d+)(\.txt)?", path.name)
        if not match:
            continue
        step = int(match[1])
        if match[2] and path.is_file():
            entry = (read_ascii_profile, path)
        elif not match[2] and (path / "particles" / "Header").is_file():
            entry = (read_native_profile, path / "particles")
        else:
            continue
        if step in found:
            raise ValueError(f"{directory}: both ASCII and native profiles for {prefix} step {step}")
        found[step] = entry
    return found


def load_profiles(directory, *, start, end, x_min, x_max, dt=None, output_interval=None,
                  prefix="sampling_fs_domain", allow_single_point=False):
    """Load eta[t, x] and the output times.

    Times come from the native output Header (exact, adaptive-dt safe). Legacy
    ASCII profiles have no time, so they need the fixed `dt` (time = step * dt)
    and, optionally, `output_interval` (steps) to check for missing files.
    Output times need not be uniform; large gaps (missing outputs) raise errors.
    """
    if not np.isfinite([start, end, x_min, x_max]).all():
        raise ValueError("Sampling parameters must be finite")
    if end <= start or x_max <= x_min:
        raise ValueError("Require increasing time/space bounds")
    if dt is not None and (not np.isfinite(dt) or dt <= 0):
        raise ValueError("dt must be positive and finite")
    if output_interval is not None and output_interval <= 0:
        raise ValueError("output_interval must be positive")
    every = []
    for step, (reader, path) in find_profiles(directory, prefix).items():
        time = read_native_time(path.parent) if reader is read_native_profile else None
        if time is None:
            if dt is None:
                raise ValueError(f"{path}: no output time recorded; pass the fixed dt")
            time = step * dt
        every.append((time, step, reader, path))
    every.sort(key=lambda item: item[0])
    all_times = np.array([item[0] for item in every])
    inside = (all_times >= start - 1e-10) & (all_times < end - 1e-10)
    snapshots = [item for item, keep in zip(every, inside) if keep]
    if len(snapshots) < 8:
        raise ValueError(f"{directory}: fewer than eight profiles in [{start}, {end}) s")
    times = np.array([item[0] for item in snapshots])
    steps = np.array([item[1] for item in snapshots])
    if np.any(np.diff(times) <= 0) or np.any(np.diff(steps) <= 0):
        raise ValueError("Duplicate or out-of-order profile times")
    if output_interval is not None and not np.all(np.diff(steps) == output_interval):
        raise ValueError("Missing/duplicate profiles or output interval does not match filenames")
    # Gap check includes the nearest outputs just outside the window.
    first, last = int(np.flatnonzero(inside)[0]), int(np.flatnonzero(inside)[-1])
    neighbourhood = all_times[max(0, first - 1):last + 2]
    spacing = np.diff(neighbourhood)
    typical = float(np.median(np.diff(times)))
    if spacing.max() > 2.5 * typical:
        gap = int(np.argmax(spacing))
        raise ValueError(f"Missing profiles: gap of {spacing[gap]:.4g} s after t={neighbourhood[gap]:.6g} s "
                         f"(median interval {typical:.4g} s)")
    covered_start = first > 0 or times[0] - start < 0.5 * typical
    covered_end = last < len(all_times) - 1 or end - times[-1] < 1.5 * typical
    if not (covered_start and covered_end):
        raise ValueError("Profiles do not cover the requested time window")

    reference = None
    elevations = []
    for _, _, reader, path in snapshots:
        count, data, _ = reader(path)
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


def resample_uniform(times, *series, rtol=1e-6):
    """Linearly interpolate (time, point) series onto a uniform time grid.

    Returns (uniform_times, resampled..., resampled_flag). The grid starts at the
    first output time with the median output interval; series already uniform
    within `rtol` are returned unchanged. Only the k-omega FFT needs this;
    accuracy norms and harmonic fits use the actual output times.
    """
    times = np.asarray(times, dtype=float)
    spacing = np.diff(times)
    step = float(np.median(spacing))
    if np.allclose(spacing, step, rtol=rtol, atol=0):
        return (times, *series, False)
    count = int(np.floor((times[-1] - times[0]) / step + 1e-9)) + 1
    grid = times[0] + step * np.arange(count)
    out = []
    for values in series:
        values = np.asarray(values, dtype=float)
        out.append(np.column_stack([np.interp(grid, times, values[:, j]) for j in range(values.shape[1])]))
    return (grid, *out, True)


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
    parser.add_argument("--dt", type=float, help="fixed dt in seconds (legacy ASCII profiles only)")
    parser.add_argument("--output-interval", type=int, help="sampler interval in steps (legacy ASCII only)")
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
    times, eta, resampled = resample_uniform(times, eta)
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
        "nt": len(times), "nx": len(x), "resampled_to_uniform_time": resampled,
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
