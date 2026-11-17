import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from temporal_analysis import (
    Record, compare, load_hos, load_sgf, make_plots, save_results,
    select_events, unique_samples, upcrossing_waves, validate_record, water_level_moments,
)


class ClockAndWaveTests(unittest.TestCase):
    def test_hos_clock_recovery_is_checked_and_uses_selected_probe(self):
        with tempfile.TemporaryDirectory() as name:
            path = Path(name) / "hos.dat"
            times = 1000 + np.arange(20) / 141.421
            text = '# f_out:: 141.421\n# Tp:: 1.6971\n# idim:: 1\nTITLE="probes"\nVARIABLES="time" "p1" "p2"\n'
            text += "".join(f"{time:.5E} {i:g} {-i:g}\n" for i, time in enumerate(times))
            path.write_text(text)
            record = load_hos(path)
            self.assertGreater(record.metadata["printed_time_ties"], 0)
            np.testing.assert_allclose(record.time, times, atol=1e-12)
            np.testing.assert_allclose(record.eta, -np.arange(20))
            path.write_text(text.replace("141.421", "100"))
            with self.assertRaisesRegex(ValueError, "rounding"):
                load_hos(path)

    def test_native_crest_is_not_reduced_by_averaging_tied_times(self):
        times = np.arange(0.1, 3.91, 0.05)
        eta = np.sin(2 * np.pi * times)
        peak = int(np.argmin(abs(times - 1.25)))
        times = np.insert(times, peak, times[peak])
        eta = np.insert(eta, peak, 2.0)
        record = Record("test", "test", times, eta, {})
        validate_record(record)
        waves = upcrossing_waves(record, 0)
        self.assertAlmostEqual(max(wave["crest_height_m"] for wave in waves), 2)
        t, z = unique_samples(record)
        self.assertAlmostEqual(np.interp(1.25, t, z), 1.5)
        self.assertEqual(len(record.time), len(times))
        self.assertAlmostEqual(record.eta.max(), 2)

    def test_complete_upcrossing_waves_exclude_partial_record_edges(self):
        times = np.linspace(0.15, 3.85, 3701)
        record = Record("test", "test", times, 0.1 * np.sin(2 * np.pi * times), {})
        waves = upcrossing_waves(record, 0)
        self.assertEqual(len(waves), 2)
        for wave, crest in zip(waves, (1.25, 2.25)):
            self.assertAlmostEqual(wave["crest_time_s"], crest)
            self.assertAlmostEqual(wave["period_s"], 1)
            self.assertAlmostEqual(wave["crest_height_m"], 0.1)
            self.assertAlmostEqual(wave["wave_height_m"], 0.2)

    def test_adaptive_sampling_does_not_bias_time_weighted_moments(self):
        times = np.array([0, 0.01, 0.02, 0.1, 0.2, 0.3, 0.5, 1.0])
        record = Record("test", "test", times, times.copy(), {})
        mean, sigma = water_level_moments(record)
        self.assertAlmostEqual(mean, 0.5)
        self.assertAlmostEqual(sigma, np.sqrt(1 / 12))

    def test_wrong_gauge_nonfinite_and_backward_timestamps_fail(self):
        with tempfile.TemporaryDirectory() as name:
            path = Path(name) / "sgf.txt"
            np.savetxt(path, np.column_stack((np.arange(8), np.ones(8), np.zeros((8, 2)))))
            with self.assertRaisesRegex(ValueError, "gauge"):
                load_sgf(path, "test", "test")
        for times, eta in ((np.r_[np.arange(7), 5], np.zeros(8)),
                           (np.arange(8), np.full(8, np.nan))):
            with self.assertRaises(ValueError):
                validate_record(Record("test", "test", times, eta, {}))

    def test_event_selection_ranks_reference_crests_and_separates_contexts(self):
        waves = [{"crest_time_s": time, "crest_height_m": height}
                 for time, height in ((1, 10), (10, 9), (11, 8), (30, 7), (40, 6))]
        events = select_events(waves, 0, 50, 3, 3)
        self.assertEqual([wave["crest_time_s"] for wave in events], [10, 30, 40])
        for wave in waves:
            wave["wave_height_m"] = wave["crest_time_s"]
        events = select_events(waves, 0, 50, 3, 3, rank_by="height")
        self.assertEqual([wave["crest_time_s"] for wave in events], [40, 30, 11])
        with self.assertRaises(ValueError):
            select_events(waves, 0, 50, 30, 3)


class ComparisonTests(unittest.TestCase):
    def example(self):
        times = np.arange(0, 20.001, 0.002)
        records = [
            Record("hos", "HOS", times, 0.1 * np.sin(2 * np.pi * times), {}),
            Record("sgf_coarse", "SGF coarse", times,
                   0.09 * np.sin(2 * np.pi * (times - 0.03)) + 0.002, {}),
            Record("sgf_fine", "SGF fine", times,
                   0.098 * np.sin(2 * np.pi * (times - 0.01)) - 0.001, {}),
        ]
        return compare(records, start=0.5, end=19.5, peak_period=1,
                       event_count=3, half_window_periods=1)

    def test_event_amplitudes_delays_and_absolute_level_bias_are_retained(self):
        report, _, _, _ = self.example()
        stats = report["statistics"]["sgf_coarse"]
        self.assertAlmostEqual(stats["mean_water_level_m"], 0.002)
        self.assertAlmostEqual(stats["maximum_crest_height_m"], 0.09)
        self.assertLess(stats["temporal_agreement"]["correlation"], 0.99)
        self.assertAlmostEqual(stats["native_maximum_elevation_m"], 0.092)
        for event in report["events"]:
            peak = event["series"]["sgf_coarse"]
            self.assertAlmostEqual(peak["crest_time_lag_s"], 0.03)
            self.assertAlmostEqual(peak["crest_height_relative_error"], -0.1)
            self.assertFalse(peak["peak_at_matching_interval_edge"])

    def test_persistent_outputs_match_statistics_and_event_shapes(self):
        report, arrays, records, waves = self.example()
        with tempfile.TemporaryDirectory() as name:
            output = Path(name)
            save_results(report, arrays, records, waves, output, plots=False)
            saved = json.loads((output / "temporal_comparison.json").read_text())
            self.assertEqual(saved, report)
            with (output / "largest_wave_events.csv").open() as stream:
                self.assertEqual(len(list(csv.DictReader(stream))), 9)
            with (output / "crest_statistics.csv").open() as stream:
                self.assertEqual(len(list(csv.DictReader(stream))), 3)
            with np.load(output / "temporal_comparison.npz") as data:
                self.assertEqual(len(data["hos_crest_height_m"]),
                                 report["statistics"]["hos"]["complete_wave_count"])
                np.testing.assert_array_equal(data["hos_eta_m"], arrays["hos_eta_m"])

    @unittest.skipUnless(importlib.util.find_spec("matplotlib"), "optional Matplotlib dependency")
    def test_event_and_crest_figures_are_saved(self):
        report, arrays, records, _ = self.example()
        with tempfile.TemporaryDirectory() as name:
            output = Path(name)
            make_plots(report, arrays, records, output)
            for name in ("largest_wave_events.png", "crest_height_statistics.png"):
                self.assertGreater((output / name).stat().st_size, 1000)


if __name__ == "__main__":
    unittest.main()
