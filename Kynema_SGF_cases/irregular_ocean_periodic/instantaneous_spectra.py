"""Compare HOS and sampled free-surface modal amplitudes in the kx-ky plane."""

from dataclasses import dataclass
from pathlib import Path
import re

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np


CASE_DIR = Path(__file__).resolve().parent
SGF_CASES = (
    ("3_ref_levels", r"SGF: $\Delta z=4.7$ m, AR=2", "sampling"),
    ("4_ref_levels", r"SGF: $\Delta z=2.3$ m, AR=2", "sampling"),
    ("4_ref_levels_AR2x", r"SGF: $\Delta z=1.2$ m, AR=4", "sampling"),
    ("5_ref_levels", r"SGF: $\Delta z=1.2$ m, AR=2", "sampling"),
    ("5_ref_levels_AR2x", r"SGF: $\Delta z=0.6$ m, AR=4", "sampling"),
    ("5_ref_levels_AR2x", "Target", "ow_sampling"),
)


@dataclass
class Spectrum:
    time: float
    kx: np.ndarray
    ky: np.ndarray
    amplitude: np.ndarray


def read_zone(path):
    with path.open() as stream:
        header = stream.readline()
    time = re.search(r'ZONE\s+T\s*=\s*"\s*([^"]+)"', header)
    nx = re.search(r"\bI\s*=\s*(\d+)", header)
    ny = re.search(r"\bJ\s*=\s*(\d+)", header)
    if time is None or nx is None or ny is None:
        raise ValueError(f"{path}: expected a Tecplot ZONE with time, I and J")
    return float(time[1]), (int(ny[1]), int(nx[1]))


def read_hos(directory):
    reference = directory / "a_3d_time0.txt"
    _, shape = read_zone(reference)
    coordinates = np.loadtxt(reference, skiprows=1, usecols=(0, 1))
    if coordinates.shape != (np.prod(shape), 2):
        raise ValueError(f"{reference}: coordinates do not match zone dimensions")
    kx = coordinates[:, 0].reshape(shape)
    ky = coordinates[:, 1].reshape(shape)
    snapshots = []
    for path in directory.glob("a_3d_time*.txt"):
        time, this_shape = read_zone(path)
        if this_shape != shape:
            raise ValueError(f"{path}: zone dimensions differ from {reference}")
        values = np.loadtxt(path, skiprows=1, ndmin=2)
        shared_coordinates = path != reference
        if shared_coordinates:
            with path.open() as stream:
                header = stream.readline()
            if re.search(r"D\s*=\s*\(\s*1\s*,\s*2\s*\)", header) is None:
                raise ValueError(f"{path}: expected shared kx/ky coordinates D=(1,2)")
        column = 0 if shared_coordinates else 2
        amplitude = values[:, column]
        if amplitude.size != np.prod(shape):
            raise ValueError(f"{path}: amplitude count does not match zone dimensions")
        if not np.all(np.isfinite(values)) or np.any(amplitude < 0):
            raise ValueError(f"{path}: expected finite, nonnegative amplitudes")
        snapshots.append(Spectrum(time, kx, ky, amplitude.reshape(shape)))
    return sorted(snapshots, key=lambda snapshot: snapshot.time)


def surface_spectrum(data, time):
    if data.ndim != 2 or data.shape[1] < 3 or not np.all(np.isfinite(data[:, :3])):
        raise ValueError("Surface samples must contain finite x, y and elevation columns")
    x, y = np.unique(data[:, 0]), np.unique(data[:, 1])
    if len(x) < 2 or len(y) < 2 or len(x) * len(y) != len(data):
        raise ValueError("Surface samples must form a complete 2D rectangular grid")
    ordered = data[np.lexsort((data[:, 0], data[:, 1]))]
    xx, yy = np.meshgrid(x, y)
    if not np.array_equal(ordered[:, :2], np.column_stack((xx.ravel(), yy.ravel()))):
        raise ValueError("Surface grid has duplicate or missing coordinate pairs")
    if not np.allclose(np.diff(x), x[1] - x[0]) or not np.allclose(
        np.diff(y), y[1] - y[0]
    ):
        raise ValueError("Surface samples must be uniformly spaced for the FFT")
    eta = ordered[:, 2].reshape(len(y), len(x))
    amplitude = np.abs(np.fft.rfft2(eta, norm="forward"))
    # HOS doubles positive-kx coefficients, except the DC and even-grid Nyquist columns.
    stop = -1 if len(x) % 2 == 0 else None
    amplitude[:, 1:stop] *= 2
    kx = 2 * np.pi * np.fft.rfftfreq(len(x), d=x[1] - x[0])
    ky = np.fft.fftshift(2 * np.pi * np.fft.fftfreq(len(y), d=y[1] - y[0]))
    return Spectrum(time, *np.meshgrid(kx, ky), np.fft.fftshift(amplitude, axes=0))


def read_surface(path):
    with path.open() as stream:
        match = re.fullmatch(r"#\s*t\s*=\s*(\S+)\s*", stream.readline())
    if match is None:
        raise ValueError(f"{path}: expected '# t = <seconds>' header")
    return surface_spectrum(np.loadtxt(path, usecols=(0, 1, 2)), float(match[1]))


def load_cases(data_dir):
    cases = [("HOS-Ocean", read_hos(data_dir / "HOS"))]
    for directory, label, prefix in SGF_CASES:
        paths = sorted((data_dir / "SGF" / directory).glob(f"{prefix}*_fs.txt"))
        if not paths:
            raise ValueError(f"No {prefix} snapshots found in {data_dir / 'SGF' / directory}")
        snapshots = sorted((read_surface(path) for path in paths), key=lambda s: s.time)
        cases.append((label, snapshots))
    count = len(cases[0][1])
    if count < 2 or any(len(snapshots) != count for _, snapshots in cases):
        raise ValueError("Comparison requires at least two snapshots and equal counts per case")
    # Match chronological snapshots without presenting different output times as identical.
    for column in zip(*(snapshots for _, snapshots in cases)):
        times = [snapshot.time for snapshot in column]
        if max(times) - min(times) > 1:
            raise ValueError(f"Snapshots are not time-aligned within 1 s: {times}")
    return cases


def make_comparison(data_dir=CASE_DIR / "data", output_dir=CASE_DIR / "plotting_outputs"):
    cases = load_cases(Path(data_dir))
    vmax = max(snapshot.amplitude.max() for _, snapshots in cases for snapshot in snapshots)
    if vmax <= 1e-6:
        raise ValueError("No amplitudes above the 1e-6 m plotting floor")
    norm = LogNorm(vmin=1e-6, vmax=vmax, clip=True)
    with plt.rc_context({"text.usetex": False}):
        fig, axes = plt.subplots(
            len(cases), len(cases[0][1]), figsize=(7, 19),
            sharex=True, sharey=True, squeeze=False, layout="constrained",
        )
        for row, (label, snapshots) in enumerate(cases):
            for col, snapshot in enumerate(snapshots):
                ax = axes[row, col]
                mesh = ax.pcolormesh(
                    snapshot.kx, snapshot.ky, snapshot.amplitude,
                    shading="auto", cmap="viridis", norm=norm, rasterized=True,
                )
                ax.set_title(f"{label}\n$t={snapshot.time:g}$ s", fontsize=10)
                ax.set_aspect("equal", adjustable="box")
                if col == 0:
                    ax.set_ylabel(r"$k_y$ [rad m$^{-1}$]")
                if row == len(cases) - 1:
                    ax.set_xlabel(r"$k_x$ [rad m$^{-1}$]")
        fig.colorbar(
            mesh, ax=axes.ravel().tolist(), label=r"$|a_\eta|$ [m]",
            extend="min", orientation="horizontal", fraction=0.015, pad=0.025,
            aspect=35,
        )
        fig.suptitle("Instantaneous 2D surface-mode amplitudes", fontsize=14)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_dir / "instantaneous_2D_spectra.png", dpi=200)
    return fig


if __name__ == "__main__":
    figure = make_comparison()
    plt.close(figure)
    print(CASE_DIR / "plotting_outputs" / "instantaneous_2D_spectra.png")
