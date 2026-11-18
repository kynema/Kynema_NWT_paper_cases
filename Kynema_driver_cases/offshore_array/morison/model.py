"""Published rigid-jacket geometry and normal-strip Morison integration."""

import hashlib
import re
from pathlib import Path

import numpy as np


REFERENCE = Path(__file__).parent / "reference/IEA-15-240-RWT-Jacket_HydroDyn.dat"
SOURCE_SHA256 = "a017b4d7711e3673e88f520dfe4f61a88e2086af10169607736319c68885549d"
ORIGINS = np.array([[-390, -390, -60], [-390, 390, -60],
                    [390, -390, -60], [390, 390, -60]], dtype=float)
CHANNELS = ("Fx", "Fy", "Fz", "Mx", "My", "Mz")


def finite(name, value, shape=None):
    value = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(value)) or (shape is not None and value.shape != shape):
        raise ValueError(f"{name}: expected finite values with shape {shape}, got {value.shape}")
    return value


def clock(time):
    time = finite("time", time)
    if time.ndim != 1 or len(time) < 2 or np.any(np.diff(time) <= 0):
        raise ValueError("time must be a strictly increasing vector with at least two samples")
    return time


def published_tables():
    """Read the pinned joint, property and member tables without changing the source."""
    data = REFERENCE.read_bytes()
    if hashlib.sha256(data).hexdigest() != SOURCE_SHA256:
        raise ValueError("Published HydroDyn source checksum does not match the pinned model")
    lines = data.decode().splitlines()

    def table(label):
        index = next(i for i, line in enumerate(lines) if re.search(r"\b" + label + r"\b", line))
        count = int(lines[index].split()[0])
        rows = []
        for line in lines[index + 1:]:
            tokens = line.split()
            if tokens and tokens[0].isdigit():
                rows.append(tokens)
            if len(rows) == count:
                return rows
        raise ValueError(f"Incomplete source table: {label}")

    return dict(joints=table("NJoints"), properties=table("NPropSets"), members=table("NMembers"))


def published_members():
    """Extract only the supported geometry tables, not a runnable HydroDyn deck."""
    tables = published_tables()
    joints = {int(row[0]): finite("joint", row[1:4], (3,)) for row in tables["joints"]}
    properties = {int(row[0]): float(row[1]) for row in tables["properties"]}
    members = []
    for row in tables["members"]:
        if row[3] != row[4] or row[6] != "1" or row[7].upper() != "FALSE":
            raise ValueError("Only untapered, simple-coefficient strip-theory members are supported")
        members.append((int(row[0]), joints[int(row[1])], joints[int(row[2])],
                        properties[int(row[3])], float(row[5])))
    return members


def build_geometry(spacing=0.5):
    if not np.isfinite(spacing) or spacing <= 0:
        raise ValueError("spacing must be positive")
    starts, ends, diameters, jackets, member_ids = [], [], [], [], []
    for jacket, origin in enumerate(ORIGINS):
        for member, a, b, diameter, source_spacing in published_members():
            a, b = a.copy(), b.copy()
            # Clip the 0.1 m buried pile ends to the CFD seabed, retaining source above it.
            if min(a[2], b[2]) < -60:
                crossing = a + (-60 - a[2]) / (b[2] - a[2]) * (b - a)
                if a[2] < -60:
                    a = crossing
                else:
                    b = crossing
            n = int(np.ceil(np.linalg.norm(b - a) / min(spacing, source_spacing)))
            points = a + np.linspace(0, 1, n + 1)[:, None] * (b - a)
            points[:, :2] += origin[:2]
            starts.extend(points[:-1])
            ends.extend(points[1:])
            diameters.extend([diameter] * n)
            jackets.extend([jacket] * n)
            member_ids.extend([member] * n)
    return dict(a=np.array(starts), b=np.array(ends), diameter=np.array(diameters),
                jacket=np.array(jackets), member=np.array(member_ids), origins=ORIGINS.copy())


def wet_geometry(geometry, eta):
    """Clip each strip against linearly interpolated endpoint surface elevations."""
    a, b = geometry["a"], geometry["b"]
    eta = finite("eta", eta, (len(a), 2))
    phi_a, phi_b = eta[:, 0] - a[:, 2], eta[:, 1] - b[:, 2]
    wet_a, wet_b = phi_a >= 0, phi_b >= 0
    crossing = np.zeros(len(a))
    mixed = wet_a != wet_b
    crossing[mixed] = phi_a[mixed] / (phi_a[mixed] - phi_b[mixed])
    lower = np.where(~wet_a & wet_b, crossing, 0)
    upper = np.where(wet_a & ~wet_b, crossing, 1)
    fraction = np.where(wet_a | wet_b, upper - lower, 0)
    xyz = a + (0.5 * (lower + upper))[:, None] * (b - a)
    length = np.linalg.norm(b - a, axis=1) * fraction
    return xyz, length


def integrate(geometry, eta, velocity, acceleration, rho=1025, cd=1.05, cm=1.2):
    """Integrate normal drag/inertia; velocity and acceleration are at wet midpoints."""
    coefficients = finite("rho, Cd, Cm", [rho, cd, cm], (3,))
    if rho <= 0 or np.any(coefficients[1:] < 0):
        raise ValueError("rho must be positive; Cd and Cm must be nonnegative")
    xyz, length = wet_geometry(geometry, eta)
    n = len(xyz)
    velocity = finite("velocity", velocity, (n, 3))
    acceleration = finite("acceleration", acceleration, (n, 3))
    tangent = geometry["b"] - geometry["a"]
    tangent /= np.linalg.norm(tangent, axis=1)[:, None]
    normal_u = velocity - np.sum(velocity * tangent, axis=1)[:, None] * tangent
    normal_a = acceleration - np.sum(acceleration * tangent, axis=1)[:, None] * tangent
    diameter = geometry["diameter"]
    drag = 0.5 * rho * cd * diameter[:, None] * np.linalg.norm(normal_u, axis=1)[:, None] * normal_u
    inertia = rho * cm * (np.pi * diameter**2 / 4)[:, None] * normal_a
    force = (drag + inertia) * length[:, None]
    arm = xyz - geometry["origins"][geometry["jacket"]]
    strips = np.column_stack((force, np.cross(arm, force)))
    loads = np.zeros((len(geometry["origins"]), 6))
    np.add.at(loads, geometry["jacket"], strips)
    return loads


def array_loads(loads, origins):
    loads = finite("loads", loads)
    origins = finite("origins", origins, (loads.shape[-2], 3))
    result = np.sum(loads, axis=-2)
    result[..., 3:] += np.sum(np.cross(origins, loads[..., :3]), axis=-2)
    return result
