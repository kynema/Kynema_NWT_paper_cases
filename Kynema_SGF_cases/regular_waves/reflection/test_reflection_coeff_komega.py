from contextlib import redirect_stderr, redirect_stdout
import io
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from reflection_coeff_komega import decompose, load_profiles, main


class DecompositionTests(unittest.TestCase):
    def waves(self, ratio=0.2, direction=1, wavelength=4.0, frequency=0.625):
        x = np.arange(256) / 32
        times = np.arange(1024) / 32
        phase = 2 * np.pi * x[None, :] / wavelength
        clock = 2 * np.pi * frequency * times[:, None]
        eta = 0.1 * np.cos(direction * phase - clock)
        eta += 0.1 * ratio * np.cos(direction * phase + clock + 0.7)
        return x, times, eta

    def test_known_reflection_and_direction(self):
        for ratio in (0, 0.01, 0.2, 1.0):
            for direction in (-1, 1):
                with self.subTest(ratio=ratio, direction=direction):
                    x, times, eta = self.waves(ratio, direction)
                    result = decompose(x, times, eta + 3, incident_direction=direction)
                    self.assertAlmostEqual(result.reflection_coefficient, ratio, places=10)
                    self.assertAlmostEqual(result.incident_variance, 0.1**2 / 2, places=10)
                    self.assertAlmostEqual(result.reflected_variance, (0.1 * ratio)**2 / 2, places=10)
                    incident = np.where(direction * result.k > 0, result.variance, 0)
                    peak = np.unravel_index(incident.argmax(), incident.shape)
                    self.assertAlmostEqual(result.k[peak[1]], direction * 2 * np.pi / 4)
                    self.assertAlmostEqual(result.omega[peak[0]], 2 * np.pi * 0.625)

    def test_harmonics_and_stationary_surface_are_excluded(self):
        x, times, eta = self.waves()
        phase = 2 * np.pi * x[None, :] / 4
        clock = 2 * np.pi * 0.625 * times[:, None]
        eta += 0.05 * np.cos(2 * phase + 2 * clock) + 0.02 * np.cos(phase)
        self.assertAlmostEqual(decompose(x, times, eta).reflection_coefficient, 0.2, places=10)

    def test_off_bin_waves_have_small_but_nonzero_leakage(self):
        x, times, eta = self.waves(0.1, wavelength=3.9, frequency=0.637)
        self.assertAlmostEqual(decompose(x, times, eta).reflection_coefficient, 0.1, delta=0.005)
        x, times, eta = self.waves(0, wavelength=3.9, frequency=0.637)
        self.assertLess(decompose(x, times, eta).reflection_coefficient, 0.01)

    def test_invalid_data_and_unresolved_bands_fail(self):
        x, times, eta = self.waves()
        irregular = times.copy()
        irregular[5] += 0.001
        nonfinite = eta.copy()
        nonfinite[0, 0] = np.nan
        for xx, tt, ee, band in (
            (x, irregular, eta, (0.4, 0.85)),
            (x, times, nonfinite, (0.4, 0.85)),
            (x, times, eta.T, (0.4, 0.85)),
            (x, times, eta, (0.626, 0.627)),
            (x, times, eta, (0, 1)),
            (x, times, np.zeros_like(eta), (0.4, 0.85)),
        ):
            with self.subTest(band=band), self.assertRaises(ValueError):
                decompose(xx, tt, ee, band)


class ProfileTests(unittest.TestCase):
    def write_profiles(self, directory):
        x = np.arange(32) / 4
        for step in range(0, 320, 10):
            t = step * 0.01
            eta = np.cos(2 * np.pi * (x / 4 - 0.625 * t))
            eta += 0.25 * np.cos(2 * np.pi * (x / 4 + 0.625 * t))
            # Sampler rows are not ordered by x; extra fields are not needed.
            rows = np.column_stack((x, np.zeros_like(x), eta, np.ones_like(x)))[::-1]
            np.savetxt(directory / f"sampling_fs_domain{step:05d}.txt", rows,
                       header="32\n0\n3\n4\n0", comments="")

    def load(self, directory):
        return load_profiles(directory, dt=0.01, output_interval=10, start=0, end=3.2,
                             x_min=0, x_max=8)

    def test_loader_and_cli_outputs(self):
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            self.write_profiles(directory)
            x, times, eta = self.load(directory)
            self.assertEqual(eta.shape, (32, 32))
            self.assertTrue(np.all(np.diff(x) > 0))
            np.testing.assert_allclose(times, np.arange(32) * 0.1)
            with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
                main([str(directory), "--dt", "0.01", "--start", "0", "--end", "3.2",
                      "--x-min", "0", "--x-max", "8",
                      "--output", str(directory / "result")])
            summary = json.loads((directory / "result.json").read_text())
            self.assertAlmostEqual(summary["reflection_coefficient"], 0.25, places=10)
            with np.load(directory / "result.npz") as saved:
                self.assertEqual(saved["variance_m2"].shape,
                                 (len(saved["omega_rad_per_s"]), len(saved["k_rad_per_m"])))

    def test_missing_snapshots_changed_grid_and_truncated_window_fail(self):
        for fault in ("missing", "grid", "start", "end", "count"):
            with self.subTest(fault=fault), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                self.write_profiles(directory)
                path = directory / "sampling_fs_domain00150.txt"
                if fault == "missing":
                    path.unlink()
                elif fault == "start":
                    (directory / "sampling_fs_domain00000.txt").unlink()
                elif fault == "end":
                    (directory / "sampling_fs_domain00310.txt").unlink()
                else:
                    text = path.read_text()
                    if fault == "count":
                        text = text.replace("32\n", "31\n", 1)
                        path.write_text(text)
                    else:
                        rows = np.loadtxt(path, skiprows=5)
                        rows[:, 0] += 0.01
                        np.savetxt(path, rows, header="32\n0\n3\n4\n0", comments="")
                with self.assertRaises(ValueError):
                    self.load(directory)


if __name__ == "__main__":
    unittest.main()
