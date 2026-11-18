#!/usr/bin/env python3
"""Export and strictly validate SeaState v5 WaveMod=6 files (13 headers, x-fast rows)."""

import argparse
import csv
import hashlib
import json
from contextlib import ExitStack
from pathlib import Path

import numpy as np

from prepare import axes


EXTENSIONS = ("Vxi", "Vyi", "Vzi", "Axi", "Ayi", "Azi", "DynP", "Elev")


def load_config(case):
    config = json.loads((Path(case) / "case.json").read_text())
    if config["openfast_version"] != "v5.0.0":
        raise ValueError("Expected the pinned OpenFAST v5.0.0 configuration")
    for name, expected in config["input_sha256"].items():
        if hashlib.sha256((Path(case) / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"{name} differs from its manifest; regenerate a consistent case")
    return config


def export(case, jacket, bundle):
    """Bundle arrays are memory-mapped .npy files, not a huge compressed NPZ."""
    config = load_config(case)
    if jacket not in range(1, 5):
        raise ValueError("jacket must be 1 through 4")
    bundle = Path(bundle)
    metadata = json.loads((bundle / "metadata.json").read_text())
    if metadata["jacket"] != jacket or metadata["registration_confirmed"] is not True:
        raise ValueError("Bundle must identify this jacket and confirm spatial/temporal registration")
    if not isinstance(metadata["source"], str) or not metadata["source"].strip():
        raise ValueError("Bundle requires a nonempty provenance source string")
    if metadata["dynamic_pressure_units"] != "Pa" or metadata["acceleration_convention"] != "eulerian":
        raise ValueError("Expected dynamic pressure in Pa and Eulerian fixed-position acceleration")
    if metadata.get("velocity_convention") != "total":
        raise ValueError("Bundle velocity_convention must be 'total': wave plus current, included exactly once")
    current = np.asarray(metadata.get("uniform_current"), dtype=float)
    if (current.shape != (3,) or not np.all(np.isfinite(current)) or
            not np.allclose(current, config["uniform_current"], rtol=0, atol=1e-12)):
        raise ValueError("Bundle uniform_current must match case.json in world-axis m/s")
    x, y, z = axes(config)
    for key, expected in (("x", x), ("y", y), ("z", z)):
        actual = np.load(bundle / f"{key}.npy", allow_pickle=False)
        if actual.shape != expected.shape or not np.allclose(actual, expected, rtol=0, atol=1e-9):
            raise ValueError(f"Bundle {key} grid does not match the generated SeaState grid")
    time = np.load(bundle / "time.npy", allow_pickle=False)
    expected_time = config["physical_time_origin"] + np.arange(config["samples"]) * config["dt"]
    if time.shape != expected_time.shape or not np.allclose(time, expected_time, rtol=0, atol=1e-8):
        raise ValueError("Bundle times do not match the declared physical clock")
    shape = (config["samples"], len(z), len(y), len(x))
    arrays = {}
    for key, expected in (("velocity", shape + (3,)), ("acceleration", shape + (3,)),
                          ("dynamic_pressure", shape), ("elevation", (shape[0], len(y), len(x)))):
        arrays[key] = np.load(bundle / f"{key}.npy", mmap_mode="r", allow_pickle=False)
        if arrays[key].shape != expected:
            raise ValueError(f"{key} must have shape {expected}, got {arrays[key].shape}")
        if arrays[key].dtype.kind not in "fiu":
            raise ValueError(f"{key} must contain real numeric data")
    output = Path(case) / "kinematics" / f"jacket{jacket}"
    if output.exists():
        raise FileExistsError(f"{output}: refusing to overwrite an existing incident realization")
    output.mkdir(parents=True)
    # Validate each frame before writing it. A failure leaves files marked incomplete.
    incomplete = output / "INCOMPLETE"
    incomplete.write_text("Export not complete; do not run HydroDyn with these files.\n")
    with ExitStack() as stack:
        streams = [stack.enter_context((output / f"waves.{ext}").open("w")) for ext in EXTENSIONS]
        for extension, stream in zip(EXTENSIONS, streams):
            header = ["SeaState v5.0.0 external incident kinematics", extension,
                      f"Source: {' '.join(metadata['source'].split())}", f"Jacket: {jacket}",
                      f"Physical time at driver t=0: {config['physical_time_origin']}",
                      f"dt: {config['dt']}", f"Samples including endpoints: {config['samples']}",
                      f"Grid nx ny nz: {len(x)} {len(y)} {len(z)}",
                      "x-fast, then y, then bottom-to-top z, then time",
                      "No time/index columns; each row contains one x sweep",
                      "Velocities m/s; accelerations m/s^2; dynamic pressure Pa; elevation m",
                      "Dynamic pressure excludes hydrostatic pressure",
                      "Data follow; exactly 13 header lines"]
            stream.write("\n".join(header) + "\n")
        for index in range(config["samples"]):
            frames = {key: value[index] for key, value in arrays.items()}
            for key, frame in frames.items():
                if not np.all(np.isfinite(frame)):
                    raise ValueError(f"Nonfinite {key} at sample {index}; export is INCOMPLETE")
            fields = ([frames["velocity"][..., i] for i in range(3)] +
                      [frames["acceleration"][..., i] for i in range(3)] +
                      [frames["dynamic_pressure"], frames["elevation"]])
            for stream, field in zip(streams, fields):
                np.savetxt(stream, field.reshape(-1, len(x)), fmt="%.9e")
    (output / "provenance.json").write_text(json.dumps(metadata, indent=2, allow_nan=False) + "\n")
    incomplete.unlink()


def validate(case):
    config = load_config(case)
    x, y, z = axes(config)
    for jacket in range(1, 5):
        root = Path(case) / "kinematics" / f"jacket{jacket}"
        if (root / "INCOMPLETE").exists():
            raise ValueError(f"{root}: incomplete export")
        for extension in EXTENSIONS:
            path = root / f"waves.{extension}"
            expected_rows = config["samples"] * len(y) * (1 if extension == "Elev" else len(z))
            with path.open() as stream:
                for _ in range(13):
                    if not stream.readline():
                        raise ValueError(f"{path}: missing required 13-line header")
                count = 0
                for count, line in enumerate(stream, 1):
                    tokens = line.partition("!")[0].split()
                    if len(tokens) != len(x):
                        raise ValueError(f"{path}: data row {count} must contain exactly {len(x)} values")
                    try:
                        values = np.array(tokens, float)
                    except ValueError:
                        raise ValueError(f"{path}: nonnumeric data at row {count}") from None
                    if not np.all(np.isfinite(values)):
                        raise ValueError(f"{path}: nonfinite data at row {count}")
                if count != expected_rows:
                    raise ValueError(f"{path}: expected {expected_rows} data rows, found {count}")
    print("All 32 external kinematics files have the required dimensions and finite numeric data.")


def grid(case, jacket, output):
    config = load_config(case)
    x, y, z = axes(config)
    origin = config["jacket_base_origins"][jacket - 1]
    with Path(output).open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["index", "x_local", "y_local", "z_SWL", "x_world", "y_world"])
        index = 0
        for zi in z:
            for yi in y:
                for xi in x:
                    writer.writerow([index, xi, yi, zi, xi + origin[0], yi + origin[1]])
                    index += 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    p = commands.add_parser("export")
    p.add_argument("case", type=Path)
    p.add_argument("--jacket", type=int, choices=range(1, 5), required=True)
    p.add_argument("--bundle", type=Path, required=True)
    p = commands.add_parser("validate")
    p.add_argument("case", type=Path)
    p = commands.add_parser("grid")
    p.add_argument("case", type=Path)
    p.add_argument("--jacket", type=int, choices=range(1, 5), required=True)
    p.add_argument("--output", type=Path, required=True)
    args = vars(parser.parse_args())
    command = args.pop("command")
    {"export": export, "validate": validate, "grid": grid}[command](**args)


if __name__ == "__main__":
    main()
