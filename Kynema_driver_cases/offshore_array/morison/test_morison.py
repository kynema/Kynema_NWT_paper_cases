import argparse
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from compare import read_cfd, welch
from hos import sample
from model import ORIGINS, array_loads, build_geometry, integrate, wet_geometry
from run import demo, geometry_id, load_geometry, prepare, run_external


def single_strip(a=(0, 0, -2), b=(0, 0, 0)):
    return dict(a=np.array([a], float), b=np.array([b], float), diameter=np.ones(1),
                jacket=np.zeros(1, int), origins=np.array([[0, 0, -2]], float))


class MorisonTests(unittest.TestCase):
    def test_partial_wetting_drag_and_moment(self):
        geometry = single_strip()
        xyz, length = wet_geometry(geometry, [[-0.5, -0.5]])
        np.testing.assert_allclose(xyz, [[0, 0, -1.25]])
        np.testing.assert_allclose(length, [1.5])
        loads = integrate(geometry, [[-0.5, -0.5]], [[2, 0, 0]], [[0, 0, 0]],
                          rho=1000, cd=1, cm=0)
        np.testing.assert_allclose(loads, [[3000, 0, 0, 0, 2250, 0]])
        reverse = single_strip(a=(0, 0, 0), b=(0, 0, -2))
        np.testing.assert_allclose(integrate(reverse, [[-0.5, -0.5]], [[2, 0, 0]],
                                             [[0, 0, 0]], rho=1000, cd=1, cm=0), loads)

    def test_inertia_coefficient_and_dry_strip(self):
        geometry = single_strip()
        loads = integrate(geometry, [[1, 1]], [[0, 0, 0]], [[1, 0, 0]], rho=1000)
        expected = 1000 * 1.2 * np.pi / 4 * 2
        np.testing.assert_allclose(loads, [[expected, 0, 0, 0, expected, 0]])
        np.testing.assert_array_equal(integrate(geometry, [[-3, -3]], [[2, 0, 0]],
                                                [[1, 0, 0]]), np.zeros((1, 6)))

    def test_inclined_member_normal_projection(self):
        geometry = single_strip(a=(0, 0, -1), b=(1, 0, 0))
        loads = integrate(geometry, [[1, 1]], [[1, 0, 1]], [[2, 0, 2]])
        np.testing.assert_allclose(loads, 0, atol=1e-10)
        loads = integrate(geometry, [[1, 1]], [[0, 0, 0]], [[1, 0, 0]], rho=1000, cm=1)
        expected = 1000 * np.pi / 4 * np.sqrt(2) * np.array([0.5, 0, -0.5])
        np.testing.assert_allclose(loads[0, :3], expected)

    def test_published_geometry_and_four_identical_transforms(self):
        geometry = build_geometry(1)
        self.assertTrue(np.all(geometry["a"][:, 2] >= -60))
        self.assertTrue(np.all(geometry["b"][:, 2] >= -60))
        self.assertTrue(np.all(geometry["diameter"][geometry["member"] == 13] == 3))
        n = len(geometry["a"])
        loads = integrate(geometry, np.zeros((n, 2)), np.tile([1, 0, 0], (n, 1)),
                          np.zeros((n, 3)))
        np.testing.assert_allclose(loads, np.tile(loads[0], (4, 1)), atol=1e-6)
        expected = loads.sum(axis=0)
        expected[3:] += np.cross(ORIGINS, loads[:, :3]).sum(axis=0)
        np.testing.assert_allclose(array_loads(loads, ORIGINS), expected)

    def test_array_moments_shift_before_summation(self):
        loads = np.array([[[0, 2, 0, 0, 0, 0], [0, 3, 0, 0, 0, 0]]], float)
        np.testing.assert_allclose(array_loads(loads, [[10, 0, 0], [20, 0, 0]]),
                                   [[0, 5, 0, 0, 0, 80]])

    def test_hos_derivative_uses_fixed_wet_midpoint_and_explicit_registration(self):
        class Field:
            end = 10

            def elevation(self, xy, time):
                return np.full(len(xy), -0.5 + 0.1 * time)

            def velocity(self, xyz, time):
                return xyz + np.array([time**2, 3 * time, 4 * time])

        geometry = single_strip()
        eta, u, a = sample(Field(), geometry, 1, [10, 20], time_shift=1, derivative_step=0.01)
        np.testing.assert_allclose(eta, [[-0.3, -0.3]])
        np.testing.assert_allclose(u, [[14, 26, 6.85]])
        np.testing.assert_allclose(a, [[4, 3, 4]], atol=1e-10)
        with self.assertRaisesRegex(ValueError, "outside"):
            sample(Field(), geometry, 0, [0, 0])

    def test_invalid_inputs_fail_explicitly(self):
        geometry = single_strip()
        for kwargs in ({"rho": 0}, {"cd": -1}, {"cm": np.nan}):
            with self.assertRaises(ValueError):
                integrate(geometry, [[0, 0]], [[1, 0, 0]], [[0, 0, 0]], **kwargs)
        with self.assertRaisesRegex(ValueError, "eta"):
            wet_geometry(geometry, [[np.nan, 0]])
        with self.assertRaises(ValueError):
            build_geometry(0)

    def test_native_force_columns_and_restarted_clock_rejection(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "forces.dat"
            path.write_text("Force output\nTime pressure viscous moment\n"
                            "0 1 2 3 10 20 30 100 200 300\n"
                            "1 2 3 4 20 30 40 200 300 400\n")
            t, loads = read_cfd(path, "nalu10", 1)
            np.testing.assert_array_equal(t, [0, 1])
            np.testing.assert_array_equal(loads[0], [11, 22, 33, 100, 200, 300])
            with path.open("a") as stream:
                stream.write("1 2 3 4 20 30 40 200 300 400\n")
            with self.assertRaisesRegex(ValueError, "strictly increasing"):
                read_cfd(path, "nalu10", 1)
            path.write_text("Force output\n0 1 broken 3 10 20 30 100 200 300\n")
            with self.assertRaisesRegex(ValueError, "invalid numeric"):
                read_cfd(path, "nalu10", 1)

    def test_cross_phase_sign_and_spectral_scaling(self):
        dt, n = 0.1, 4096
        time = np.arange(n) * dt
        frequency = 10 / (512 * dt)
        a = np.cos(2 * np.pi * frequency * time)[:, None]
        b = 2 * np.cos(2 * np.pi * frequency * time + np.pi / 6)[:, None]
        result = welch(a, b, dt)
        index = np.argmin(abs(result["frequency"] - frequency))
        self.assertAlmostEqual(result["phase_deg"][index, 0], 30, places=2)
        self.assertAlmostEqual(result["cfd_psd"][index, 0] / result["reference_psd"][index, 0], 4, places=4)
        self.assertGreater(result["coherence"][index, 0], 0.999)
        self.assertAlmostEqual(np.trapezoid(result["reference_psd"][:, 0], result["frequency"]), 0.5, places=4)

    def test_external_kinematics_and_full_comparison_smoke_test(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            geometry = prepare(output / "prepared", 1)
            np.testing.assert_array_equal(load_geometry(output / "prepared/geometry.npz")["a"], geometry["a"])
            n = len(geometry["a"])
            kinematics = output / "kinematics.npz"
            np.savez_compressed(kinematics, time=[0, 1], eta=np.zeros((2, n, 2)),
                                velocity=np.tile([1, 0, 0], (2, n, 1)),
                                acceleration=np.zeros((2, n, 3)), source="analytic constant current",
                                acceleration_convention="eulerian", geometry_sha256=geometry_id(geometry))
            args = argparse.Namespace(geometry=output / "prepared/geometry.npz", kinematics=kinematics,
                                      output=output / "loads.npz", cd=1.05, cm=1.2)
            run_external(args)
            with np.load(args.output) as data:
                np.testing.assert_allclose(data["loads"][0], data["loads"][1])
            demo(argparse.Namespace(output=output / "demo", plot=False))
            report = json.loads((output / "demo/comparison/comparison.json").read_text())
            self.assertEqual(len(report["statistics"]), 30)
            self.assertTrue(report["largest_load_events"])
            for row in report["statistics"]:
                if row["model_std"] > 1e-5:
                    self.assertAlmostEqual(row["cfd_std"] / row["model_std"], 1.1, places=8)
                    self.assertAlmostEqual(row["demeaned_rmse"] / row["model_std"], 0.1, places=8)
                    self.assertAlmostEqual(row["correlation"], 1, places=8)


if __name__ == "__main__":
    unittest.main()
