import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np

from collect import collect, read_output
from prepare import axes, generate
from wavekin import export, validate
from compare import compare
from model import CHANNELS, published_members


def make_bundle(path, config, jacket, varying=False):
    path.mkdir()
    x, y, z = axes(config)
    time = config["physical_time_origin"] + np.arange(config["samples"]) * config["dt"]
    shape = (len(time), len(z), len(y), len(x))
    velocity = np.zeros(shape + (3,))
    velocity[..., 0] = config["uniform_current"][0] * jacket
    if varying:
        velocity[..., 0] = np.arange(np.prod(shape)).reshape(shape)
    for name, value in dict(x=x, y=y, z=z, time=time, velocity=velocity,
                            acceleration=np.zeros_like(velocity),
                            dynamic_pressure=np.zeros(shape),
                            elevation=np.zeros((len(time), len(y), len(x)))).items():
        np.save(path / f"{name}.npy", value)
    (path / "metadata.json").write_text(json.dumps(dict(
        jacket=jacket, registration_confirmed=True, source="SYNTHETIC test field, not demonstration data",
        uniform_current=config["uniform_current"], velocity_convention="total",
        dynamic_pressure_units="Pa", acceleration_convention="eulerian")))


def run_driver(executable, case, name):
    result = subprocess.run([executable, str(case / name)], capture_output=True, text=True, timeout=60)
    if result.returncode or "HydroDyn Driver terminated normally." not in result.stdout:
        raise AssertionError(result.stdout + result.stderr)


class HydroDynTests(unittest.TestCase):
    def test_external_layout_and_inclusive_clock(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            case, bundle = root / "case", root / "bundle"
            config = generate(case, duration=0.3, dt=0.15, physical_time_origin=240,
                              nx_half=2, ny_half=2, nz=3)
            make_bundle(bundle, config, 1, varying=True)
            metadata = json.loads((bundle / "metadata.json").read_text())
            metadata["source"] = "SYNTHETIC source\nwith a multiline description"
            (bundle / "metadata.json").write_text(json.dumps(metadata))
            export(case, 1, bundle)
            path = case / "kinematics/jacket1/waves.Vxi"
            data = np.loadtxt(path, skiprows=13)
            expected = np.arange(3 * 3 * 3 * 3).reshape(-1, 3)
            np.testing.assert_array_equal(data, expected)
            with self.assertRaises(FileExistsError):
                export(case, 1, bundle)
            with self.assertRaisesRegex(FileExistsError, "already contains"):
                generate(case, duration=0.3, dt=0.15, nx_half=2, ny_half=2, nz=3)
            self.assertEqual(config["samples"], 3)
            self.assertIn("3 NSteps", (case / "jacket1.dvr").read_text())

    def test_bad_registration_clock_and_nonfinite_data_fail(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            case, bundle = root / "case", root / "bundle"
            config = generate(case, duration=0.3, dt=0.15, nx_half=2, ny_half=2, nz=3)
            make_bundle(bundle, config, 1)
            np.save(bundle / "time.npy", np.array([0, 0.15, 0.3]))
            with self.assertRaisesRegex(ValueError, "physical clock"):
                export(case, 1, bundle)
            np.save(bundle / "time.npy", np.array([240, 240.15, 240.3]))
            a = np.load(bundle / "velocity.npy")
            a[0, 0, 0, 0, 0] = np.nan
            np.save(bundle / "velocity.npy", a)
            with self.assertRaisesRegex(ValueError, "Nonfinite"):
                export(case, 1, bundle)
            self.assertTrue((case / "kinematics/jacket1/INCOMPLETE").exists())

    def test_all_files_validate_and_extra_rows_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            case = root / "case"
            config = generate(case, duration=0.3, dt=0.15, nx_half=2, ny_half=2, nz=3)
            for jacket in range(1, 5):
                bundle = root / f"bundle{jacket}"
                make_bundle(bundle, config, jacket)
                export(case, jacket, bundle)
            validate(case)
            with (case / "kinematics/jacket1/waves.Elev").open("a") as stream:
                stream.write("0 0 0\n")
            with self.assertRaisesRegex(ValueError, "data rows"):
                validate(case)

    def test_time_and_grid_validation(self):
        with tempfile.TemporaryDirectory() as directory:
            for kwargs in ({"duration": 1, "dt": 0.3}, {"nz": 1}, {"dt": 0}, {"member_ds": -1},
                           {"current_x": np.nan}):
                with self.assertRaises(ValueError):
                    generate(directory, **kwargs)

    def test_current_contract_rejects_wave_only_or_mismatched_bundles(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            case, bundle = root / "case", root / "bundle"
            config = generate(case, duration=0.3, dt=0.15, nx_half=2, ny_half=2, nz=3,
                              current_x=-1.5)
            make_bundle(bundle, config, 1)
            metadata = json.loads((bundle / "metadata.json").read_text())
            for change in ({"velocity_convention": "wave-only"}, {"velocity_convention": None},
                           {"uniform_current": [0, 0, 0]}, {"uniform_current": [1.5, 0, 0]},
                           {"uniform_current": [-1.5]}, {"uniform_current": None},
                           {"uniform_current": [np.nan, 0, 0]}):
                with self.subTest(change=change):
                    (bundle / "metadata.json").write_text(json.dumps({**metadata, **change}))
                    with self.assertRaisesRegex(ValueError, "velocity_convention|uniform_current"):
                        export(case, 1, bundle)
                    self.assertFalse((case / "kinematics").exists())
            (bundle / "metadata.json").write_text(json.dumps(metadata))
            export(case, 1, bundle)
            values = np.loadtxt(case / "kinematics/jacket1/waves.Vxi", skiprows=13)
            np.testing.assert_allclose(values, -1.5)

    @unittest.skipUnless(os.environ.get("HYDRODYN_DRIVER"), "Set HYDRODYN_DRIVER for executable smoke tests")
    def test_actual_driver_regular_and_external_four_jackets(self):
        executable = os.environ["HYDRODYN_DRIVER"]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            case = root / "case"
            config = generate(case, duration=2.4, dt=0.15, physical_time_origin=240,
                              nx_half=3, ny_half=3, nz=9)
            run_driver(executable, case, "smoke_regular.dvr")
            time, loads = read_output(case / "smoke_regular.out")
            self.assertEqual(len(time), 246)
            self.assertAlmostEqual(time[-1], 36.75)
            self.assertGreater(loads[:, 0].std(), 0)
            summary = (case / "smoke_regular.HD.sum").read_text().splitlines()
            index = next(i for i, line in enumerate(summary)
                         if line.split()[:1] == ["MemberID"] and "FilledFlag" in line)
            header = summary[index].split()
            rows = [line.split() for line in summary[index + 2:index + 59]]
            flooded = [int(row[0]) for row in rows if row[header.index("FilledFlag")] == "T"]
            self.assertEqual(flooded, list(range(1, 16)))
            radii = {member: diameter / 2 for member, _, _, diameter, _ in published_members()}
            for row in rows:
                self.assertAlmostEqual(float(row[header.index("R1")]), radii[int(row[0])])
            for jacket in range(1, 5):
                bundle = root / f"bundle{jacket}"
                make_bundle(bundle, config, jacket)
                export(case, jacket, bundle)
                run_driver(executable, case, f"jacket{jacket}.dvr")
            validate(case)
            collect(case, root / "loads.npz")
            with np.load(root / "loads.npz", allow_pickle=False) as data:
                np.testing.assert_allclose(data["time"], 240 + np.arange(17) * 0.15)
                fx = data["loads"][0, :, 0]
                self.assertTrue(np.all(fx > 0))
                self.assertTrue(np.all(data["loads"][0, :, 4] > 0))
                np.testing.assert_allclose(fx / fx[0], np.arange(1, 5)**2, rtol=1e-6)
                time, raw = read_output(case / "jacket1.out")
                np.testing.assert_allclose(data["loads"][:, 0, 4], raw[:, 4] + 60 * raw[:, 0])
                np.testing.assert_allclose(data["loads"][:, 0, 3], raw[:, 3] - 60 * raw[:, 1])
                paths = {}
                for jacket in range(4):
                    paths[jacket + 1] = root / f"synthetic_cfd{jacket + 1}.csv"
                    np.savetxt(paths[jacket + 1], np.column_stack((data["time"], data["loads"][:, jacket])),
                               delimiter=",", header="time," + ",".join(CHANNELS), comments="")
            report = compare(root / "loads.npz", paths, root / "comparison", "csv", 1, start=240, dt=0.15)
            self.assertNotIn("Normal-strip loads only", report["warning"])
            self.assertTrue(all(row["rmse"] < 1e-6 for row in report["statistics"]))


if __name__ == "__main__":
    unittest.main()
