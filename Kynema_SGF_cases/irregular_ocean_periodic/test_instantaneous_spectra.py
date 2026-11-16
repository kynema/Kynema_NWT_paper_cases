from pathlib import Path
import tempfile
import unittest

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from instantaneous_spectra import CASE_DIR, load_cases, make_comparison, read_hos, surface_spectrum


class SurfaceSpectrumTests(unittest.TestCase):
    def samples(self, nx=16, ny=12):
        x, y = np.meshgrid(np.arange(nx) * 0.5, np.arange(ny) * 0.75)
        eta = 3 * np.cos(2 * np.pi * (2 * x / (nx * 0.5) - 3 * y / (ny * 0.75)))
        return np.column_stack((x.ravel(), y.ravel(), eta.ravel()))

    def test_oblique_wave_amplitude_and_signed_wavenumbers(self):
        for nx in (15, 16):
            with self.subTest(nx=nx):
                samples = self.samples(nx=nx)
                spectrum = surface_spectrum(samples[::-1], 12.5)
                peak = np.unravel_index(spectrum.amplitude.argmax(), spectrum.amplitude.shape)
                self.assertAlmostEqual(spectrum.amplitude[peak], 3)
                self.assertAlmostEqual(spectrum.kx[peak], 2 * np.pi * 2 / (nx * 0.5))
                self.assertAlmostEqual(spectrum.ky[peak], -2 * np.pi * 3 / (12 * 0.75))
                self.assertEqual(spectrum.time, 12.5)
                self.assertEqual(spectrum.amplitude.shape, (12, nx // 2 + 1))

    def test_dc_and_nyquist_are_not_doubled(self):
        samples = self.samples()
        samples[:, 2] = 2 + 4 * np.cos(np.pi * samples[:, 0] / 0.5)
        spectrum = surface_spectrum(samples, 0)
        zero_y = np.argmin(abs(spectrum.ky[:, 0]))
        self.assertAlmostEqual(spectrum.amplitude[zero_y, 0], 2)
        self.assertAlmostEqual(spectrum.amplitude[zero_y, -1], 4)

    def test_zero_kx_transverse_wave_is_not_doubled(self):
        samples = self.samples()
        samples[:, 2] = 6 * np.cos(2 * np.pi * 2 * samples[:, 1] / (12 * 0.75))
        spectrum = surface_spectrum(samples, 0)
        self.assertAlmostEqual(spectrum.amplitude.max(), 3)

    def test_odd_grid_last_kx_is_doubled(self):
        samples = self.samples(nx=15)
        samples[:, 2] = 4 * np.cos(2 * np.pi * 7 * samples[:, 0] / (15 * 0.5))
        spectrum = surface_spectrum(samples, 0)
        self.assertAlmostEqual(spectrum.amplitude[:, -1].max(), 4)

    def test_invalid_grids_are_rejected(self):
        missing = self.samples()[:-1]
        duplicate = self.samples()
        duplicate[1, :2] = duplicate[0, :2]
        irregular = self.samples()
        irregular[irregular[:, 0] == 0, 0] = 0.1
        nonfinite = self.samples()
        nonfinite[0, 2] = np.nan
        for samples in (missing, duplicate, irregular, nonfinite):
            with self.subTest(samples=samples.shape):
                with self.assertRaises(ValueError):
                    surface_spectrum(samples, 0)


class HOSReaderTests(unittest.TestCase):
    def test_single_zone_and_shared_coordinates(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            initial = path / "a_3d_time0.txt"
            initial.write_text('ZONE T = "0", I=2, J=2\n0 -1 1\n1 -1 2\n0 0 3\n1 0 4\n')
            self.assertEqual(len(read_hos(path)), 1)
            later = path / "a_3d_time1000.txt"
            later.write_text('ZONE T = "999.6", D=(1,2), I=2, J=2\n5\n6\n7\n8\n')
            snapshots = read_hos(path)
            self.assertEqual([s.time for s in snapshots], [0, 999.6])
            np.testing.assert_array_equal(snapshots[0].kx, snapshots[1].kx)
            np.testing.assert_array_equal(snapshots[1].amplitude, [[5, 6], [7, 8]])
            later.write_text('ZONE T = "999.6", D=(1,2), I=2, J=2\n5\n6\n7\n')
            with self.assertRaisesRegex(ValueError, "amplitude count"):
                read_hos(path)


class RepositoryComparisonTests(unittest.TestCase):
    def test_all_cases_times_and_shared_plot_scales(self):
        cases = load_cases(CASE_DIR / "data")
        self.assertEqual(len(cases), 7)
        self.assertTrue(all(len(snapshots) == 2 for _, snapshots in cases))
        self.assertEqual([s.time for s in cases[0][1]], [0, 999.6])
        self.assertEqual([s.time for s in cases[-1][1]], [0, 999.989])
        self.assertTrue(all(s.amplitude.shape == (256, 129) for _, ss in cases for s in ss))
        with tempfile.TemporaryDirectory() as directory:
            figure = make_comparison(output_dir=directory)
            try:
                self.assertEqual(len(figure.axes), 15)
                panels = figure.axes[:-1]
                norm = panels[0].collections[0].norm
                self.assertEqual(norm.vmin, 1e-6)
                for ax in panels:
                    self.assertIs(ax.collections[0].norm, norm)
                    self.assertEqual(ax.get_xlim(), panels[0].get_xlim())
                    self.assertEqual(ax.get_ylim(), panels[0].get_ylim())
                    self.assertEqual(ax.get_aspect(), 1)
                self.assertIn("999.6", panels[1].get_title())
                self.assertIn("999.989", panels[-1].get_title())
                self.assertGreater((Path(directory) / "instantaneous_2D_spectra.png").stat().st_size, 0)
            finally:
                plt.close(figure)


if __name__ == "__main__":
    unittest.main()
