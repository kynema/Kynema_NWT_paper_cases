#!/usr/bin/env python3
"""Convert standalone HydroDyn outputs to base-referenced SI loads for comparison."""

import argparse
import json
from pathlib import Path

import numpy as np

from prepare import ORIGINS, SOURCE_SHA256
from wavekin import load_config


HYDRO_CHANNELS = ("HydroFxi", "HydroFyi", "HydroFzi", "HydroMxi", "HydroMyi", "HydroMzi")


def read_output(path):
    with Path(path).open() as stream:
        for line in stream:
            names = line.split()
            if names and names[0] == "Time":
                break
        else:
            raise ValueError(f"{path}: no Time/channel header found")
        units = stream.readline().split()
        if len(units) != len(names) or "INVALID" in units:
            raise ValueError(f"{path}: invalid or missing channel units")
        indices = [names.index(name) for name in HYDRO_CHANNELS]
        if units[0] != "(s)" or [units[i] for i in indices] != ["(N)"] * 3 + ["(N-m)"] * 3:
            raise ValueError(f"{path}: expected seconds, N and N-m")
        data = np.loadtxt(stream, ndmin=2)
    if data.shape[1] != len(names) or not np.all(np.isfinite(data)):
        raise ValueError(f"{path}: nonfinite or incomplete output data")
    if len(data) < 2 or np.any(np.diff(data[:, 0]) <= 0):
        raise ValueError(f"{path}: output times must be strictly increasing")
    return data[:, 0], data[:, indices]


def collect(case, output):
    case = Path(case)
    config = load_config(case)
    histories = [read_output(case / f"jacket{jacket}.out") for jacket in range(1, 5)]
    expected_time = np.arange(config["samples"]) * config["dt"]
    loads = []
    for time, data in histories:
        if time.shape != expected_time.shape or not np.allclose(time, expected_time, rtol=0, atol=5.1e-5):
            raise ValueError("HydroDyn times do not match case.json; check completion and driver dt")
        base_loads = data.copy()
        # HydroM channels are at local SWL origin, not the driver's PtfmRefzt.
        base_loads[:, 3:] += np.cross([0, 0, 60], base_loads[:, :3])
        loads.append(base_loads)
    metadata = dict(source="HydroDyn v5.0.0", reference_label="HydroDyn",
                    model="standalone fixed-jacket HydroDyn Morison with hydrostatics and axial loads",
                    source_sha256=SOURCE_SHA256, density=1025, Cd=1.05, Ca=0.2, Cp=1,
                    flooding=config["flooding"], marine_growth=False, wave_stretching=config["wave_stretching"],
                    uniform_current=config["uniform_current"], velocity_convention="total",
                    input_sha256=config["input_sha256"],
                    physical_time_origin=config["physical_time_origin"],
                    exclusions=["structural response", "diffraction", "inter-jacket shielding",
                                "exact CAD-joint intersection volumes", "resolved interior flooded flow"],
                    scope_warning="HydroDyn flooded-leg, bare-member model: verify CFD flooding, end conditions, "
                                  "joint overlaps, acceleration convention and vertical stretching.",
                    case_directory=str(case.resolve()))
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wb") as stream:
        np.savez_compressed(stream, time=expected_time + config["physical_time_origin"],
                            loads=np.stack(loads, axis=1), origins=ORIGINS,
                            metadata=json.dumps(metadata, allow_nan=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    collect(args.case, args.output)


if __name__ == "__main__":
    main()
