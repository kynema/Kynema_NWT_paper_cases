from contextlib import redirect_stdout
import io
import json
import math
from pathlib import Path
import tempfile
import unittest

import numpy as np

from analyze_case import accuracy_metrics, analyze, main
from generate_cases import case_specs, generate, render_input, render_refinement


def parameters(text):
    return {
        key.strip(): value.split("#")[0].strip()
        for line in text.splitlines() if "=" in line and not line.startswith("#")
        for key, value in [line.split("=", 1)]
    }


class SweepTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cases = case_specs()

    def test_sweeps_share_one_baseline_and_cover_requested_values(self):
        self.assertEqual(len({case["id"] for case in self.cases}), 24)
        for family in ("linear", "stokes5"):
            cases = [case for case in self.cases if case["family"] == family]
            mesh = [case for case in cases if "mesh" in case["sweeps"]]
            aspect = [case for case in cases if "aspect_ratio" in case["sweeps"]]
            cfl = [case for case in cases if "cfl" in case["sweeps"]]
            self.assertEqual([case["mesh_target_cells_per_height"] for case in mesh], [1, 2.5, 5, 10, 20])
            self.assertEqual([case["aspect_ratio"] for case in aspect], [1, 2, 4, 8])
            self.assertEqual(sorted(case["nominal_advective_cfl"] for case in cfl), [0.1, 0.2, 0.5, 0.8, 1])
            shared = [case for case in cases if len(case["sweeps"]) == 3]
            self.assertEqual(len(shared), 1)
            self.assertEqual(shared[0]["mesh_target_cells_per_height"], 5)
            self.assertEqual(len({case["dx_m"] for case in aspect + cfl}), 1)
            self.assertEqual(len({case["dy_m"] for case in aspect + cfl}), 1)

    def test_rendered_geometry_resolution_refinement_and_sampling(self):
        for case in self.cases:
            with self.subTest(case=case["id"]):
                p = parameters(render_input(case))
                lo, hi = (np.fromstring(p[key], sep=" ") for key in ("geometry.prob_lo", "geometry.prob_hi"))
                n = np.fromstring(p["amr.n_cell"], sep=" ")
                spacing = (hi - lo) / n / 2**int(p["amr.max_level"])
                self.assertAlmostEqual(spacing[0], spacing[1])
                self.assertAlmostEqual(spacing[0] / spacing[2], case["aspect_ratio"])
                cph = float(p["OceanWaves.Wave1.wave_height"]) / spacing[2]
                expected = case["mesh_target_cells_per_height"] * case["aspect_ratio"]
                self.assertLess(abs(cph / expected - 1), 0.05)
                wavelength = float(p["OceanWaves.Wave1.wave_length"])
                self.assertEqual(p["zlo.type"], '"slip_wall"')
                self.assertEqual(p["zhi.type"], '"slip_wall"')
                self.assertEqual(float(p["OceanWaves.Wave1.numerical_beach_length"]), wavelength)
                self.assertEqual(float(p["OceanWaves.Wave1.numerical_beach_length_factor"]), 2)
                self.assertEqual(p["OceanWaves.Wave1.initialize_wave_field"], "true")
                self.assertEqual(float(p["OceanWaves.Wave1.timeramp_period"]), 0)
                self.assertGreaterEqual((hi[0] - lo[0]) / wavelength, 10)
                self.assertAlmostEqual(-lo[2], float(p["OceanWaves.Wave1.water_depth"]))
                self.assertGreaterEqual(case["analysis_time_s"][1] - case["analysis_time_s"][0], 40 * 1.59)
                self.assertGreater(float(p["time.stop_time"]), case["analysis_time_s"][1])
                self.assertNotIn("time.initial_dt", p)
                self.assertAlmostEqual(float(p["time.fixed_dt"]), case["dt_s"])
                for group, label, phase in (
                    ("sampling_fs_domain", "fs_domain", "vof"),
                    ("sampling_fs_target", "fs_target", "ow_vof"),
                ):
                    self.assertEqual(p.get(f"{group}.{label}.phase_field", "vof"), phase)
                    self.assertEqual(int(p.get(f"{group}.output_frequency", "10")),
                                     case["output_interval_steps"])
                    self.assertLess(case["sample_dt_s"], 1 / (2 * case["frequency_band_hz"][1]))
                boxes = render_refinement(case).splitlines()
                for row in (2, 4):
                    box = np.fromstring(boxes[row], sep=" ")
                    self.assertLess(box[2], -case["wave_height_m"])
                    self.assertGreater(box[5], case["wave_height_m"])
                    self.assertLessEqual(box[2], -0.5 * case["wave_height_m"])
                    self.assertTrue(np.all(box[:3] >= lo))
                    self.assertTrue(np.all(box[3:] <= hi))

    def test_default_controls_are_omitted_but_physics_remains_explicit(self):
        omitted = (
            "time.max_step", "time.plot_interval", "time.checkpoint_interval",
            "incflo.use_godunov", "incflo.godunov_type", "incflo.diffusione_type",
            "MultiPhase.verbose", "amr.ref_ratio", "amr.blocking_factor_x",
            "amr.regrid_int", "sampling_fs_domain.fs_domain.phase_field",
        )
        for case in self.cases:
            with self.subTest(case=case["id"]):
                p = parameters(render_input(case))
                for key in omitted:
                    self.assertNotIn(key, p)
                self.assertEqual(p["incflo.diffusion_type"], "0")
                self.assertEqual(float(p.get("time.cfl", "1")), case["nominal_advective_cfl"])
                self.assertEqual("time.cfl" in p, case["nominal_advective_cfl"] != 1)
                for group in ("sampling_fs_domain", "sampling_fs_target"):
                    self.assertEqual(f"{group}.output_frequency" in p,
                                     case["output_interval_steps"] != 10)
                self.assertEqual(p["incflo.gravity"], "0 0 -9.81")
                self.assertEqual(float(p["transport.viscosity_fluid1"]), 0)
                self.assertEqual(float(p["transport.viscosity_fluid2"]), 0)
                self.assertEqual(float(p["OceanWaves.Wave1.zero_sea_level"]), 0)
                self.assertEqual(float(p["OceanWaves.Wave1.wave_phase_offset_radians"]), 0)

    def test_timesteps_preserve_cfl_and_match_solver_vof_convention(self):
        phase = np.linspace(0, 2 * np.pi, 65536, endpoint=False)
        for family in ("linear", "stokes5"):
            cases = [case for case in self.cases if case["family"] == family]
            mesh = [case for case in cases if "mesh" in case["sweeps"]]
            scaling = [case["dt_s"] / case["dx_m"] for case in mesh]
            np.testing.assert_allclose(scaling, scaling[0], rtol=1e-12)
            cfl = [case for case in cases if "cfl" in case["sweeps"]]
            scaling = [case["dt_s"] / case["nominal_advective_cfl"] for case in cfl]
            np.testing.assert_allclose(scaling, scaling[0], rtol=1e-12)
            aspect = [case for case in cases if "aspect_ratio" in case["sweeps"]]
            self.assertTrue(np.all(np.diff([case["dt_s"] for case in aspect]) < 0))
            for case in cases:
                self.assertAlmostEqual(case["dt_s"] * case["cfl_rate_estimate_per_s"],
                                       case["nominal_advective_cfl"])
                if family == "linear":
                    a = case["wave_height_m"] / 2
                    k = 2 * np.pi / case["wave_length_m"]
                    d = case["water_depth_m"]
                    omega = np.sqrt(9.81 * k * np.tanh(k * d))
                    eta = a * np.cos(phase)
                    u = a * omega * np.cosh(k * (d + eta)) / np.sinh(k * d) * np.cos(phase)
                    w = a * omega * np.sinh(k * (d + eta)) / np.sinh(k * d) * np.sin(phase)
                    expected = 2 * np.max(abs(u) / case["dx_m"] + abs(w) / case["dz_m"])
                    self.assertAlmostEqual(case["dt_s"] * expected,
                                           case["nominal_advective_cfl"], delta=1e-6)

    def test_generation_is_reproducible_and_linear_steepness_is_bounded(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            generate(root)
            first = {path.relative_to(root): path.read_bytes() for path in root.rglob("*") if path.is_file()}
            generate(root)
            second = {path.relative_to(root): path.read_bytes() for path in root.rglob("*") if path.is_file()}
            self.assertEqual(first, second)
            self.assertEqual(len(list(root.rglob("input.inp"))), 24)
            for path in root.rglob("case.json"):
                case = json.loads(path.read_text())
                self.assertTrue((path.parent / "static_box.txt").is_file())
                if case["family"] == "linear":
                    self.assertAlmostEqual(case["ka"], 0.05)
                    self.assertGreater(case["wave_height_m"], 0.0256)


class AnalysisTests(unittest.TestCase):
    def write_run(self, root):
        case = next(case for case in case_specs() if case["id"] == "linear/mesh_cph5")
        case.update(dt_s=0.01, output_interval_steps=10, sample_dt_s=0.1,
                    analysis_time_s=[0.0, 3.2], sample_points=57)
        (root / "case.json").write_text(json.dumps(case))
        post = root / "post_processing"
        post.mkdir()
        x = np.linspace(8, 36, 57)
        a = case["wave_height_m"] / 2
        for step in range(0, 320, 10):
            phase = 2 * np.pi * x / 4
            clock = 2 * np.pi * 0.625 * step * 0.01
            target = a * np.cos(phase - clock)
            eta = target + 0.2 * a * np.cos(phase + clock + 0.7)
            for prefix, z in (("sampling_fs_domain", eta), ("sampling_fs_target", target)):
                rows = np.column_stack((x, np.zeros_like(x), z))[::-1]
                np.savetxt(post / f"{prefix}{step:05d}.txt", rows,
                           header="57\n0\n3\n4\n0", comments="")
        return case

    def test_unified_analysis_and_persistent_outputs(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            self.write_run(root)
            summary, arrays = analyze(root)
            self.assertAlmostEqual(summary["reflection_coefficient"], 0.2, places=10)
            self.assertAlmostEqual(summary["accuracy"]["l2_over_height"],
                                   0.2 / (2 * math.sqrt(2)), places=10)
            self.assertLess(summary["target_apparent_reflection_coefficient"], 1e-10)
            self.assertAlmostEqual(summary["harmonic_analysis"]["fundamental_frequency_hz"], 0.625, places=7)
            self.assertEqual(len(summary["harmonic_analysis"]["harmonics"]), 5)
            self.assertEqual(arrays["harmonic_phase_error_rad"].shape, (5, 56))
            self.assertAlmostEqual(np.mean(arrays["rmse_by_x_m"]**2),
                                   summary["accuracy"]["l2_m"]**2)
            with redirect_stdout(io.StringIO()):
                main([str(root)])
            self.assertEqual(json.loads((root / "analysis/results.json").read_text()), summary)
            with np.load(root / "analysis/results.npz") as saved:
                np.testing.assert_allclose(saved["variance_m2"], arrays["variance_m2"])

    def test_missing_target_is_an_error_not_an_accuracy_fallback(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            self.write_run(root)
            (root / "post_processing/sampling_fs_target00150.txt").unlink()
            with self.assertRaisesRegex(ValueError, "Missing"):
                analyze(root)

    def test_phase_error_is_not_aligned_away_and_invalid_inputs_fail(self):
        phase = np.linspace(0, 2 * np.pi, 100, endpoint=False)
        target = np.cos(phase)
        result = accuracy_metrics(np.cos(phase - np.pi / 2), target, 2)
        self.assertAlmostEqual(result["l2_m"], 1)
        for eta, ref, height in ((target[:-1], target, 2), (target * np.nan, target, 2),
                                 (target, target, 0)):
            with self.assertRaises(ValueError):
                accuracy_metrics(eta, ref, height)


if __name__ == "__main__":
    unittest.main()
