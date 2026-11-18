#!/usr/bin/env python3
"""Prepare, run and compare a published-diameter, fixed-jacket Morison benchmark."""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from compare import compare
from model import CHANNELS, SOURCE_SHA256, build_geometry, clock, finite, integrate, wet_geometry


def geometry_id(geometry):
    digest = hashlib.sha256()
    for key in ("a", "b", "diameter", "jacket", "member", "origins"):
        dtype = "<i8" if key in ("jacket", "member") else "<f8"
        digest.update(np.asarray(geometry[key], dtype=dtype).tobytes())
    return digest.hexdigest()


def load_geometry(path):
    with np.load(path, allow_pickle=False) as data:
        geometry = {key: data[key].copy() for key in
                    ("a", "b", "diameter", "jacket", "member", "origins")}
        expected = str(data["geometry_sha256"])
    if expected != geometry_id(geometry):
        raise ValueError("Geometry checksum mismatch; regenerate using prepare")
    return geometry


def prepare(output, spacing):
    geometry = build_geometry(spacing)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    identity = geometry_id(geometry)
    np.savez_compressed(output / "geometry.npz", **geometry, geometry_sha256=identity)
    with (output / "strip_endpoints.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["strip", "jacket", "member", "endpoint", "x", "y", "z"])
        for i in range(len(geometry["a"])):
            for endpoint, key in enumerate(("a", "b")):
                writer.writerow([i, int(geometry["jacket"][i]) + 1,
                                 int(geometry["member"][i]), endpoint, *geometry[key][i]])
    metadata = dict(geometry="published diameters; no diameter scaling or marine growth",
                    source_sha256=SOURCE_SHA256, geometry_sha256=identity,
                    strip_count=len(geometry["a"]), max_spacing=spacing,
                    world_datum="SWL z=0; seabed z=-60; no jacket rotation",
                    warning="Existing CFD mesh has approximately doubled tube diameters. "
                            "Rebuild and audit it before treating a comparison as validation.")
    (output / "geometry.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return geometry


def save_loads(path, geometry, time, loads, source, cd, cm, **details):
    time = clock(time)
    loads = finite("loads", loads, (len(time), 4, 6))
    metadata = dict(source=source, source_sha256=SOURCE_SHA256,
                    geometry_sha256=geometry_id(geometry), geometry="published diameters",
                    model="fixed normal-strip Morison; midpoint quadrature; no growth",
                    acceleration_convention="Eulerian partial du/dt at fixed world coordinates",
                    density=1025, Cd=cd, Cm=cm,
                    exclusions=["hydrostatics", "axial/end loads", "joint-overlap corrections",
                                "structural dynamics", "diffraction", "inter-jacket shielding"],
                    **details)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, time=time, loads=loads, origins=geometry["origins"],
                        metadata=json.dumps(metadata, allow_nan=False))


def run_external(args):
    geometry = load_geometry(args.geometry)
    with np.load(args.kinematics, allow_pickle=False) as data:
        if str(data["geometry_sha256"]) != geometry_id(geometry):
            raise ValueError("Kinematics do not reference the prepared geometry")
        if str(data["acceleration_convention"]) != "eulerian":
            raise ValueError("This benchmark requires acceleration_convention='eulerian'")
        source = str(data["source"])
        if not source.strip() or source == "synthetic":
            raise ValueError("External kinematics require a nonempty, nonsynthetic provenance string")
        time = clock(data["time"])
        n = len(geometry["a"])
        eta = finite("eta", data["eta"], (len(time), n, 2))
        velocity = finite("velocity", data["velocity"], (len(time), n, 3))
        acceleration = finite("acceleration", data["acceleration"], (len(time), n, 3))
    loads = np.array([integrate(geometry, eta[i], velocity[i], acceleration[i], cd=args.cd, cm=args.cm)
                      for i in range(len(time))])
    save_loads(args.output, geometry, time, loads, source, args.cd, args.cm,
               kinematics_file=str(Path(args.kinematics).resolve()))


def run_hos(args):
    from hos import Grid2Grid, sample

    if not args.registration_confirmed:
        raise ValueError("Verify HOS/CFD position, clock, direction and current, then pass --registration-confirmed")
    if not all(np.isfinite(v) for v in (args.start, args.stop, args.dt)) or args.dt <= 0:
        raise ValueError("start/stop must be finite and dt positive")
    if args.stop <= args.start:
        raise ValueError("stop must be after start")
    current_x = float(finite("current_x", args.current_x, ()))
    geometry = load_geometry(args.geometry)
    time = clock(args.start + np.arange(int(np.floor((args.stop - args.start) / args.dt)) + 1) * args.dt)
    provider = Grid2Grid(args.library, args.modes, args.zmax, args.nz)
    loads = []
    for t in time:
        eta, velocity, acceleration = sample(provider, geometry, t, args.hos_xy_shift,
                                             args.hos_time_shift, args.derivative_step,
                                             current_x=current_x)
        loads.append(integrate(geometry, eta, velocity, acceleration, cd=args.cd, cm=args.cm))
    save_loads(args.output, geometry, time, loads, "HOS via Grid2Grid", args.cd, args.cm,
               modes_file=str(Path(args.modes).resolve()), library=str(Path(args.library).resolve()),
               hos_xy_shift=args.hos_xy_shift, hos_time_shift=args.hos_time_shift,
               derivative_step=args.derivative_step, zmax=args.zmax, nz=args.nz,
               uniform_current=[current_x, 0, 0], velocity_convention="total",
               registration_confirmed=True)


def synthetic_loads(geometry, time):
    omega, depth, amplitude = 2 * np.pi / 12.25, 60, 1
    lower, upper = 0.00001, 1
    for _ in range(60):
        k = (lower + upper) / 2
        if 9.81 * k * np.tanh(k * depth) > omega**2:
            upper = k
        else:
            lower = k
    direction = np.array([np.cos(0.3), np.sin(0.3)])
    endpoints = np.stack((geometry["a"], geometry["b"]), axis=1)
    records = []
    for t in time:
        eta = amplitude * np.cos(k * np.sum(endpoints[..., :2] * direction, axis=-1) - omega * t)
        xyz, _ = wet_geometry(geometry, eta)
        theta = k * np.sum(xyz[:, :2] * direction, axis=1) - omega * t
        ch = np.cosh(k * (xyz[:, 2] + depth)) / np.sinh(k * depth)
        sh = np.sinh(k * (xyz[:, 2] + depth)) / np.sinh(k * depth)
        u = np.column_stack((amplitude * omega * (ch * np.cos(theta))[:, None] * direction,
                             amplitude * omega * sh * np.sin(theta)))
        a = np.column_stack((amplitude * omega**2 * (ch * np.sin(theta))[:, None] * direction,
                             -amplitude * omega**2 * sh * np.cos(theta)))
        records.append(integrate(geometry, eta, u, a))
    return np.array(records)


def demo(args):
    output = Path(args.output)
    geometry = prepare(output / "prepared", 1)
    time = np.arange(0, 120.01, 0.5)
    loads = synthetic_loads(geometry, time)
    save_loads(output / "synthetic_morison.npz", geometry, time, loads, "synthetic", 1.05, 1.2,
               warning="Analytic linear-wave smoke test; CFD inputs below are fabricated, not simulation data")
    paths = {}
    for jacket in range(4):
        paths[jacket + 1] = output / f"synthetic_cfd_jacket{jacket + 1}.csv"
        np.savetxt(paths[jacket + 1], np.column_stack((time, loads[:, jacket] * 1.1)),
                   delimiter=",", header="time," + ",".join(CHANNELS), comments="")
    compare(output / "synthetic_morison.npz", paths, output / "comparison",
            "csv", 1, start=0, dt=0.5, plot=args.plot)
    print("SYNTHETIC ONLY: fabricated CFD loads = 1.1 * Morison loads; no real CFD comparison.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    p = commands.add_parser("prepare", help="Export published geometry and strip endpoint requests")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--spacing", type=float, default=0.5)
    p = commands.add_parser("integrate", help="Integrate externally sampled incident kinematics")
    p.add_argument("geometry", type=Path)
    p.add_argument("kinematics", type=Path)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--cd", type=float, default=1.05)
    p.add_argument("--cm", type=float, default=1.2)
    p = commands.add_parser("hos", help="Sample exact HOS realization using a compatible Grid2Grid library")
    p.add_argument("geometry", type=Path)
    p.add_argument("--library", type=Path, required=True)
    p.add_argument("--modes", type=Path, required=True)
    p.add_argument("--hos-xy-shift", type=float, nargs=2, required=True,
                   help="HOS xy = CFD world xy + these explicitly verified shifts [m]")
    p.add_argument("--hos-time-shift", type=float, default=0)
    p.add_argument("--current-x", type=float, required=True,
                   help="Uniform x-current [m/s] added once to wave-only HOS velocities; demo: 1.5")
    p.add_argument("--registration-confirmed", action="store_true")
    p.add_argument("--start", type=float, default=240)
    p.add_argument("--stop", type=float, required=True)
    p.add_argument("--dt", type=float, default=0.15)
    p.add_argument("--derivative-step", type=float, default=0.02)
    p.add_argument("--zmax", type=float, default=16)
    p.add_argument("--nz", type=int, default=100)
    p.add_argument("--cd", type=float, default=1.05)
    p.add_argument("--cm", type=float, default=1.2)
    p.add_argument("--output", type=Path, required=True)
    p = commands.add_parser("compare", help="Compare four CFD histories with model loads")
    p.add_argument("loads", type=Path)
    p.add_argument("--cfd", action="append", nargs=2, required=True, metavar=("JACKET", "FILE"))
    p.add_argument("--cfd-format", choices=("csv", "nalu10"), required=True)
    p.add_argument("--force-on-structure-sign", type=int, choices=(-1, 1), required=True)
    p.add_argument("--start", type=float, default=240)
    p.add_argument("--dt", type=float, default=0.15)
    p.add_argument("--plot", action="store_true")
    p.add_argument("--output", type=Path, required=True)
    p = commands.add_parser("check-mesh", help="Audit wall tube radii in an Exodus mesh")
    p.add_argument("mesh", type=Path)
    p.add_argument("--output", type=Path, required=True)
    p = commands.add_parser("demo", help="Run a synthetic, explicitly fabricated end-to-end smoke test")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    if args.command == "prepare":
        prepare(args.output, args.spacing)
    elif args.command == "integrate":
        run_external(args)
    elif args.command == "hos":
        run_hos(args)
    elif args.command == "compare":
        ids = [int(jacket) for jacket, _ in args.cfd]
        if len(set(ids)) != len(ids):
            parser.error("Duplicate --cfd jacket IDs; stitch restart histories explicitly first")
        compare(args.loads, {int(j): Path(path) for j, path in args.cfd}, args.output,
                args.cfd_format, args.force_on_structure_sign, args.start, args.dt, args.plot)
    elif args.command == "check-mesh":
        from mesh import audit
        report = audit(args.mesh)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        print(json.dumps(report, indent=2))
        if report["status"] == "mismatch":
            raise SystemExit(2)
    elif args.command == "demo":
        demo(args)


if __name__ == "__main__":
    main()
