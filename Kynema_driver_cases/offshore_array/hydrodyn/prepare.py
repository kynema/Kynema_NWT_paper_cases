#!/usr/bin/env python3
"""Generate fixed, flooded-leg jacket decks for standalone OpenFAST v5.0.0."""

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "morison"))
from model import ORIGINS, SOURCE_SHA256, published_members, published_tables


def render(template, values):
    text = (HERE / "templates" / template).read_text()
    for key, value in values.items():
        text = text.replace("@" + key + "@", str(value))
    if re.search(r"@[A-Z_]+@", text):
        raise ValueError(f"Unresolved template token in {template}")
    return text


def axes(config):
    x = np.linspace(-25, 25, 2 * config["nx_half"] - 1)
    y = np.linspace(-25, 25, 2 * config["ny_half"] - 1)
    theta = np.arange(config["nz"] - 1, -1, -1) * np.pi / (2 * (config["nz"] - 1))
    z = -(1 - np.cos(theta)) * 60
    z[0], z[-1] = -60, 0
    return x, y, z


def generate(output, duration=3435, dt=0.15, physical_time_origin=240,
             nx_half=6, ny_half=6, nz=31, member_ds=0.5, current_x=1.5):
    if not np.all(np.isfinite([duration, dt, physical_time_origin, member_ds])):
        raise ValueError("Times and member spacing must be finite")
    if duration <= 0 or dt <= 0 or physical_time_origin < 0 or member_ds <= 0:
        raise ValueError("Duration, dt and member spacing must be positive; time origin nonnegative")
    if not np.isfinite(current_x):
        raise ValueError("Uniform x-current must be finite")
    if any(not isinstance(n, int) or n < 2 for n in (nx_half, ny_half, nz)):
        raise ValueError("Grid counts must be integers at least two")
    intervals = round(duration / dt)
    if intervals < 1 or not np.isclose(intervals * dt, duration, rtol=0, atol=1e-8):
        raise ValueError("Duration must be an integer multiple of dt")
    tables = published_tables()
    published_members()
    joints = []
    for row in tables["joints"]:
        row = row.copy()
        row[3] = f"{max(float(row[3]), -60):.10g}"
        joints.append(" ".join(row))
    members = []
    for row in tables["members"]:
        members.append(" ".join(row[:5] + ["1", "0", f"{min(member_ds, float(row[5])):.10g}",
                                          row[6], "1", row[7]]))
    output = Path(output)
    if (output / "kinematics").exists():
        raise FileExistsError("Case already contains incident data; generate a new directory to avoid stale mappings")
    output.mkdir(parents=True, exist_ok=True)
    (output / "LICENSE").write_bytes((HERE.parent / "morison/reference/LICENSE").read_bytes())
    (output / "HydroDyn.dat").write_text(render("HydroDyn.dat", dict(
        JOINTS="\n".join(joints), PROPERTIES="\n".join(" ".join(r) for r in tables["properties"]),
        MEMBERS="\n".join(members))))
    config = dict(openfast_version="v5.0.0", source_sha256=SOURCE_SHA256,
                  geometry_source="mmrocze2/IEA-15-240-RWT@da77a393569b6c5ddf7aa424e5b0a92d37eba833",
                  format_reference="OpenFAST/r-test@dd5feaaaa500ba7283140107806300d551cff0a7",
                  duration=intervals * dt, dt=dt, samples=intervals + 1,
                  physical_time_origin=physical_time_origin, nx_half=nx_half,
                  ny_half=ny_half, nz=nz, member_ds=member_ds,
                  jacket_base_origins=ORIGINS.tolist(), flooding="members 1-15; FillFSLoc=0; FillDens=1025",
                  marine_growth=False, wave_stretching="vertical; unchanged below SWL, constant above SWL",
                  uniform_current=[current_x, 0, 0],
                  hydrodyn_moment_origin=[0, 0, 0], exported_moment_origin="each jacket base at z=-60")
    x, y, z = axes(config)
    config["grid"] = dict(x=x.tolist(), y=y.tolist(), z=z.tolist(), order="time,z,y,x; x fastest")
    config["wave_arrays_minimum_bytes_per_jacket"] = (intervals + 1) * len(x) * len(y) * (7 * len(z) + 1) * 4
    for jacket in range(1, 5):
        title = f"Jacket {jacket}: same incident realization sampled at its world location; external data required."
        (output / f"jacket{jacket}_SeaState.dat").write_text(render("SeaState.dat", dict(
            TITLE=title, NX=nx_half, NY=ny_half, NZ=nz, WAVEMOD=6, DURATION=f"{intervals * dt:.12g}",
            CURRENT_DESCRIPTION=f"INCLUDED IN EXTERNAL VELOCITIES: Ux={current_x:+g} m/s; DO NOT ADD TWICE",
            DT=f"{dt:.12g}", KINROOT=f"kinematics/jacket{jacket}/waves")))
        (output / f"jacket{jacket}.dvr").write_text(render("driver.inp", dict(
            TITLE=title, SEASTATE=f"jacket{jacket}_SeaState.dat", ROOT=f"jacket{jacket}",
            NSAMPLES=intervals + 1, DT=f"{dt:.12g}")))
    (output / "smoke_SeaState.dat").write_text(render("SeaState.dat", dict(
        TITLE="SYNTHETIC regular-wave parser/run smoke test; not demonstration data.",
        CURRENT_DESCRIPTION="NONE: WAVE-ONLY SYNTHETIC SMOKE TEST",
        NX=3, NY=3, NZ=9, WAVEMOD=1, DURATION=36.75, DT=0.15, KINROOT="unused")))
    (output / "smoke_regular.dvr").write_text(render("driver.inp", dict(
        TITLE="SYNTHETIC regular-wave smoke test; fixed published-diameter flooded-leg jacket.",
        SEASTATE="smoke_SeaState.dat", ROOT="smoke_regular", NSAMPLES=246, DT=0.15)))
    config["input_sha256"] = {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(output.iterdir()) if path.suffix in (".dat", ".dvr")
    }
    (output / "case.json").write_text(json.dumps(config, indent=2, allow_nan=False) + "\n")
    return config


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "input_files")
    parser.add_argument("--duration", type=float, default=3435)
    parser.add_argument("--dt", type=float, default=0.15)
    parser.add_argument("--physical-time-origin", type=float, default=240)
    parser.add_argument("--nx-half", type=int, default=6)
    parser.add_argument("--ny-half", type=int, default=6)
    parser.add_argument("--nz", type=int, default=31)
    parser.add_argument("--member-ds", type=float, default=0.5)
    parser.add_argument("--current-x", type=float, default=1.5,
                        help="Expected uniform x-current in external total velocities [m/s]; demo: 1.5")
    args = parser.parse_args()
    config = generate(**vars(args))
    print(f"Generated four external-wave cases for {config['openfast_version']}; no real kinematics supplied.")
    print(f"Minimum seven-field+elevation array storage per jacket: "
          f"{config['wave_arrays_minimum_bytes_per_jacket'] / 1024**3:.2f} GiB, excluding other allocations.")


if __name__ == "__main__":
    main()
