import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from harmonic_analysis import harmonic_analysis, plot_harmonics
from reflection.reflection_coeff_komega import load_profiles


class HarmonicTests(unittest.TestCase):
    def samples(self, lag=0.023, ratio=0.9):
        times = 123.4 + np.arange(1000) * 0.02
        x = np.arange(8) * 0.5
        frequency = 0.637
        amplitudes = np.array([0.2, 0.025, 0.004, 0.001, 0.0004])
        target = np.full((len(times), len(x)), 0.3)
        eta = target.copy()
        for h, amplitude in enumerate(amplitudes, 1):
            phase = h * 2 * np.pi * x[None, :] / 4 + 0.17 * h
            clock = h * 2 * np.pi * frequency * (times - times[0])[:, None]
            target += amplitude * np.cos(phase - clock)
            eta += ratio * amplitude * np.cos(phase - clock + h * 2 * np.pi * frequency * lag)
        return x, times, eta, target, frequency, amplitudes

    def test_joint_off_bin_harmonics_recover_amplitude_phase_and_lag(self):
        x, times, eta, target, frequency, amplitudes = self.samples()
        summary, arrays = harmonic_analysis(times, eta, target)
        self.assertAlmostEqual(summary["fundamental_frequency_hz"], frequency, places=8)
        self.assertEqual(arrays["harmonic_phase_valid"].shape, (5, len(x)))
        for index, row in enumerate(summary["harmonics"]):
            with self.subTest(order=row["order"]):
                self.assertAlmostEqual(row["target_amplitude_rms_m"], amplitudes[index], places=9)
                self.assertAlmostEqual(row["relative_amplitude_error"], -0.1, places=7)
                self.assertAlmostEqual(row["equivalent_time_lag_s"], 0.023, places=7)
                self.assertEqual(row["phase_status"], "resolved")
        self.assertLess(summary["signal_fit_residual_rms_m"], 1e-8)
        self.assertLess(summary["target_fit_residual_rms_m"], 1e-8)

    def test_wrapped_phase_and_circular_averaging_do_not_cancel_at_180(self):
        times = np.arange(640) * 0.025
        clock = 2 * np.pi * 0.625 * times[:, None]
        target = np.cos(-clock + np.zeros((1, 2)))
        eta = np.cos(-clock + np.deg2rad([[179, -179]]))
        summary, arrays = harmonic_analysis(times, eta, target, count=1, frequency=0.625)
        self.assertGreater(abs(summary["harmonics"][0]["phase_error_circular_mean_deg"]), 179)
        np.testing.assert_allclose(np.rad2deg(arrays["harmonic_phase_error_rad"]), [[179, -179]])
        reference = np.cos(-clock + np.deg2rad([[179]]))
        shifted = np.cos(-clock + np.deg2rad([[-179]]))
        summary, _ = harmonic_analysis(times, shifted, reference, count=1, frequency=0.625)
        self.assertAlmostEqual(summary["harmonics"][0]["phase_error_circular_mean_deg"], 2)

    def test_weak_harmonics_and_vanishing_signal_have_no_reported_phase(self):
        times = np.arange(400) * 0.02
        phase = 2 * np.pi * 0.625 * times[:, None]
        target = np.cos(phase)
        eta = 0.8 * target + 0.05 * np.cos(2 * phase + 1)
        summary, arrays = harmonic_analysis(times, eta, target, frequency=0.625)
        self.assertAlmostEqual(summary["harmonics"][1]["signal_amplitude_rms_m"], 0.05)
        for row in summary["harmonics"][1:]:
            self.assertIsNone(row["phase_error_circular_mean_deg"])
            self.assertIsNone(row["relative_amplitude_error"])
            self.assertEqual(row["phase_status"], "weak_reference_or_signal")
        self.assertTrue(np.isnan(arrays["harmonic_phase_error_rad"][1:]).all())
        json.dumps(summary, allow_nan=False)
        summary, _ = harmonic_analysis(times, np.zeros_like(target), target, frequency=0.625)
        self.assertIsNone(summary["harmonics"][0]["equivalent_time_lag_s"])

    def test_signal_frequency_error_is_not_fitted_away(self):
        _, times, _, target, frequency, _ = self.samples()
        clock = 2 * np.pi * (frequency + 0.002) * (times - times[0])[:, None]
        eta = 0.2 * np.cos(np.arange(8)[None, :] * np.pi / 4 + 0.17 - clock)
        summary, _ = harmonic_analysis(times, eta, target)
        self.assertAlmostEqual(summary["fundamental_frequency_hz"], frequency, places=8)
        self.assertGreater(abs(summary["harmonics"][0]["phase_error_circular_mean_deg"]), 5)
        self.assertGreater(summary["signal_fit_residual_rms_m"], 1e-3)

    def test_invalid_fits_fail_explicitly(self):
        _, times, eta, target, _, _ = self.samples()
        for kwargs in ({"count": 1000}, {"count": 50}, {"phase_floor_ratio": 0},
                       {"frequency": np.nan}, {"frequency_band": (0.6, 0.6)}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                harmonic_analysis(times, eta, target, **kwargs)
        with self.assertRaises(ValueError):
            harmonic_analysis(times, eta, np.zeros_like(target))
        with self.assertRaises(ValueError):
            harmonic_analysis(times, eta[:-1], target)
        eta[0, 0] = np.nan
        with self.assertRaises(ValueError):
            harmonic_analysis(times, eta, target)

    def test_single_point_loader_requires_explicit_opt_in(self):
        with tempfile.TemporaryDirectory() as name:
            path = Path(name)
            for step in range(0, 160, 10):
                (path / f"sampling_fs_gauge{step:05d}.txt").write_text(
                    f"1\n0\n3\n4\n0\n10 0 {np.cos(step * 0.01)}\n"
                )
            options = dict(dt=0.001, output_interval=10, start=0, end=0.16,
                           x_min=9.99, x_max=10.01, prefix="sampling_fs_gauge")
            with self.assertRaisesRegex(ValueError, "single-point"):
                load_profiles(path, **options)
            x, times, eta = load_profiles(path, allow_single_point=True, **options)
            self.assertEqual(eta.shape, (16, 1))
            np.testing.assert_allclose(x, [10])
            np.testing.assert_allclose(times, np.arange(16) * 0.01)

    @unittest.skipUnless(importlib.util.find_spec("matplotlib"), "optional Matplotlib dependency")
    def test_harmonic_plot_is_saved_with_masked_phase(self):
        x, times, eta, target, _, _ = self.samples()
        summary, arrays = harmonic_analysis(times, eta, target)
        arrays["x_m"] = x
        arrays["harmonic_phase_valid"][1:] = False
        arrays["harmonic_phase_error_rad"][1:] = np.nan
        for row in summary["harmonics"][1:]:
            row["phase_error_circular_mean_deg"] = None
        with tempfile.TemporaryDirectory() as name:
            output = Path(name) / "harmonics.png"
            plot_harmonics({"case": "synthetic", "harmonic_analysis": summary}, arrays, output)
            self.assertGreater(output.stat().st_size, 1000)


if __name__ == "__main__":
    unittest.main()
