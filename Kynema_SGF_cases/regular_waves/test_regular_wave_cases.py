from contextlib import redirect_stdout
import io
import json
import math
from pathlib import Path
import tempfile
import unittest

import numpy as np

from analyze_case import accuracy_metrics, analyze, main
from generate_cases import (
    AIR_DENSITY, AIR_VISCOSITY, CHECKPOINT_TIME_S, FAMILIES, WATER_DENSITY, WATER_VISCOSITY,
    case_specs, generate, render_input, render_refinement,
)


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
                # Adaptive time stepping: CFL-limited dt, calibrated initial estimate
                self.assertEqual(float(p["time.fixed_dt"]), -1)
                self.assertAlmostEqual(float(p["time.initial_dt"]), case["initial_dt_s"])
                self.assertEqual(case["time_stepping"], "adaptive")
                for group, label, phase in (
                    ("sampling_fs_domain", "fs_domain", "vof"),
                    ("sampling_fs_target", "fs_target", "ow_vof"),
                ):
                    self.assertEqual(p.get(f"{group}.{label}.phase_field", "vof"), phase)
                    self.assertNotIn(f"{group}.output_frequency", p)
                    self.assertNotIn(f"{group}.output_interval", p)
                    self.assertEqual(float(p[f"{group}.output_time_interval"]),
                                     case["sample_time_interval_s"])
                    # Fifth harmonic of the band edge stays below Nyquist
                    self.assertLess(case["sample_time_interval_s"],
                                    1 / (2 * 5 * case["frequency_band_hz"][1]))
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
                self.assertNotIn("incflo.diffusion_type", p)
                # Infrequent restart checkpoints (~CHECKPOINT_TIME_S of simulated time).
                self.assertEqual(float(p["time.checkpoint_time_interval"]), CHECKPOINT_TIME_S)
                self.assertEqual(case["checkpoint_time_interval_s"], CHECKPOINT_TIME_S)
                # Native AMG bottom solver: the default bicgstab bottom solve stalls the nodal
                # projection on the two-cell periodic span (linear/mesh_cph1 hit 200 iterations).
                self.assertEqual(p["nodal_proj.bottom_solver"], "amg")
                self.assertEqual(p["mac_proj.bottom_solver"], "amg")
                # time.cfl is the adaptive limit, so it is always explicit
                self.assertEqual(float(p["time.cfl"]), case["nominal_advective_cfl"])
                for group in ("sampling_fs_domain", "sampling_fs_target"):
                    self.assertNotIn(f"{group}.output_frequency", p)
                self.assertEqual(p["incflo.gravity"], "0 0 -9.81")
                self.assertEqual(float(p["transport.viscosity_fluid1"]), WATER_VISCOSITY)
                self.assertEqual(float(p["transport.viscosity_fluid2"]), AIR_VISCOSITY)
                self.assertEqual(float(p["MultiPhase.density_fluid1"]), WATER_DENSITY)
                self.assertEqual(float(p["MultiPhase.density_fluid2"]), AIR_DENSITY)
                # wave_generation ghost fill takes the gas density from xlo.density
                self.assertEqual(float(p["xlo.density"]), AIR_DENSITY)
                self.assertEqual(float(p["OceanWaves.Wave1.zero_sea_level"]), 0)
                self.assertEqual(float(p["OceanWaves.Wave1.wave_phase_offset_radians"]), 0)

    def test_timesteps_preserve_cfl_and_match_solver_vof_convention(self):
        phase = np.linspace(0, 2 * np.pi, 65536, endpoint=False)
        for family in ("linear", "stokes5"):
            cases = [case for case in self.cases if case["family"] == family]
            mesh = [case for case in cases if "mesh" in case["sweeps"]]
            # With a fixed well-resolved factor, dt scales with dx apart from the
            # stability multipliers of poorly resolved meshes.
            scaling = [case["initial_dt_s"] * case["stability_multiplier"] / case["dx_m"] for case in mesh]
            np.testing.assert_allclose(scaling, scaling[0], rtol=1e-12)
            for case in mesh:
                self.assertEqual(case["stability_multiplier"],
                                 FAMILIES[family]["stability_multipliers"].get(
                                     case["mesh_target_cells_per_height"], 1.0))
                if case["mesh_target_cells_per_height"] >= 10:
                    self.assertEqual(case["stability_multiplier"], 1.0)
            cfl = [case for case in cases if "cfl" in case["sweeps"]]
            scaling = [case["initial_dt_s"] / case["nominal_advective_cfl"] for case in cfl]
            np.testing.assert_allclose(scaling, scaling[0], rtol=1e-12)
            aspect = [case for case in cases if "aspect_ratio" in case["sweeps"]]
            self.assertTrue(np.all(np.diff([case["initial_dt_s"] for case in aspect]) < 0))
            for case in cases:
                factor = FAMILIES[family]["cfl_rate_factor"] * case["stability_multiplier"]
                self.assertEqual(case["cfl_rate_factor"], FAMILIES[family]["cfl_rate_factor"])
                self.assertAlmostEqual(case["effective_cfl_rate_factor"], factor)
                self.assertAlmostEqual(
                    case["initial_dt_s"] * factor * case["cfl_rate_estimate_per_s"],
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
                    self.assertAlmostEqual(case["initial_dt_s"] * factor * expected,
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
        # Legacy fixed-dt ASCII run: time = step * dt
        case.update(time_stepping="fixed", dt_s=0.01, output_interval_steps=10, sample_dt_s=0.1,
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

    @staticmethod
    def write_native_output(post, prefix, step, time, x, z):
        """One AMReX native sampler output with its time in the directory Header."""
        out = post / f"{prefix}{step:05d}"
        level = out / "particles" / "Level_0"
        level.mkdir(parents=True)
        (out / "Header").write_text(f"HyperCLaw-V1.1\n1\nvof\n3\n{float(time)!r}\n")
        order = np.arange(len(x))[::-1]
        ints = np.zeros((len(x), 5), dtype="<i4")
        ints[:, 4] = order
        reals = np.column_stack((x[order], np.zeros(len(x)), z[order], np.ones(len(x)))).astype("<f8")
        (level / "DATA_00000").write_bytes(ints.tobytes() + reals.tobytes())
        (out / "particles" / "Header").write_text("\n".join(
            ["Version_Two_Dot_Zero_double", "3", "1", "vof", "3", "uid", "set_id", "probe_id",
             "1", str(len(x)), str(len(x) + 1), "0", "1", f"0 {len(x)} 0"]) + "\n")

    def test_adaptive_nonuniform_native_outputs(self):
        # Adaptive dt with time-interval output: irregular output times read from Headers.
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            case = next(case for case in case_specs() if case["id"] == "linear/mesh_cph5")
            case.update(analysis_time_s=[0.0, 32.0], sample_points=57)
            (root / "case.json").write_text(json.dumps(case))
            post = root / "post_processing"
            post.mkdir()
            x = np.linspace(8, 36, 57)
            a = case["wave_height_m"] / 2
            n = np.arange(700)
            times = n * 0.05 + 0.012 * np.sin(1.7 * n) * (n > 0)
            for step, t in zip(n * 7, times):
                phase = 2 * np.pi * x / 4
                clock = 2 * np.pi * 0.625 * t
                target = a * np.cos(phase - clock)
                eta = target + 0.2 * a * np.cos(phase + clock + 0.7)
                self.write_native_output(post, "sampling_fs_domain", step, t, x, eta)
                self.write_native_output(post, "sampling_fs_target", step, t, x, target)
            summary, _ = analyze(root)
            self.assertTrue(summary["output_times"]["resampled_to_uniform_for_komega"])
            self.assertGreater(summary["output_times"]["max_interval_s"],
                               summary["output_times"]["min_interval_s"] + 0.01)
            self.assertAlmostEqual(summary["reflection_coefficient"], 0.2, delta=2e-3)
            self.assertLess(summary["target_apparent_reflection_coefficient"], 2e-3)
            self.assertAlmostEqual(summary["harmonic_analysis"]["fundamental_frequency_hz"], 0.625, places=6)
            self.assertAlmostEqual(summary["accuracy"]["l2_over_height"],
                                   0.2 / (2 * math.sqrt(2)), delta=5e-4)

    def test_stuck_sampler_points_are_an_error(self):
        # A FreeSurfaceSampler point owned by no cell keeps a constant fallback height
        # (domain bottom or a coarse-cell face) in both the surface and the target.
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            self.write_run(root)
            for path in (root / "post_processing").glob("sampling_fs_*.txt"):
                rows = np.loadtxt(path, skiprows=5)
                rows[rows[:, 0] == 20.0, 2] = -4.0
                np.savetxt(path, rows, header="57\n0\n3\n4\n0", comments="")
            with self.assertRaisesRegex(ValueError, "never change"):
                analyze(root)

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
