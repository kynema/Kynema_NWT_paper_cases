"""Optional Exodus wall-node radius audit; not a full CAD equivalence check."""

import hashlib

import numpy as np

from model import SOURCE_SHA256, published_members


def audit(path):
    from netCDF4 import Dataset, chartostring

    side_nodes = {
        5: [[0, 1, 4], [1, 2, 4], [2, 3, 4], [0, 4, 3], [0, 3, 2, 1]],
        4: [[0, 1, 3], [1, 2, 3], [0, 3, 2], [0, 2, 1]],
    }
    with Dataset(path) as mesh:
        names = list(chartostring(mesh["ss_names"][:]))
        wall_index = names.index("wall") + 1
        xyz = np.array([mesh[key][:] for key in ("coordx", "coordy", "coordz")]).T
        elements = mesh[f"elem_ss{wall_index}"][:] - 1
        sides = mesh[f"side_ss{wall_index}"][:] - 1
        node_ids, offset = [], 0
        for block in range(1, len(mesh.dimensions["num_el_blk"]) + 1):
            connection = mesh[f"connect{block}"][:] - 1
            kind = mesh[f"connect{block}"].elem_type.upper()
            if kind not in ("TETRA", "TETRA4", "PYRAMID", "PYRAMID5"):
                raise ValueError(f"Unsupported Exodus element type {kind}")
            selected = (elements >= offset) & (elements < offset + len(connection))
            for side, local in enumerate(side_nodes[connection.shape[1]]):
                pick = selected & (sides == side)
                node_ids.extend(connection[elements[pick] - offset][:, local].ravel())
            offset += len(connection)
    wall = xyz[np.unique(node_ids)]
    if not len(wall):
        raise ValueError("No wall nodes found")
    sample = wall[::max(1, len(wall) // 10000)]
    members = published_members()
    a = np.array([row[1] for row in members]) + [0, 0, 60]
    b = np.array([row[2] for row in members]) + [0, 0, 60]
    diameter = np.array([row[3] for row in members])
    tangent = b - a
    length = np.linalg.norm(tangent, axis=1)
    tangent /= length[:, None]
    relative = sample[:, None, :] - a
    along = np.einsum("nmi,mi->nm", relative, tangent)
    radial = np.linalg.norm(relative - along[:, :, None] * tangent, axis=2)
    interior = (along > 0.5) & (along < length - 0.5)
    residuals = {}
    for scale in (1, 2):
        distance = np.min(np.where(interior, abs(radial - scale * diameter / 2), np.inf), axis=1)
        distance = distance[np.isfinite(distance)]
        if not len(distance):
            raise ValueError("No mesh points project into published member interiors")
        residuals[str(scale)] = dict(zip(("p10", "median", "p90", "p95"),
                                         np.quantile(distance, [0.1, 0.5, 0.9, 0.95]).tolist()))
    with path.open("rb") as stream:
        checksum = hashlib.sha256()
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            checksum.update(chunk)
    return dict(mesh_name=path.name, mesh_sha256=checksum.hexdigest(),
                source_sha256=SOURCE_SHA256, wall_nodes=len(wall),
                wall_bounds=[wall.min(axis=0).tolist(), wall.max(axis=0).tolist()],
                source_to_mesh_local_shift=[0, 0, 60],
                surface_distance_m_by_diameter_scale=residuals,
                status="mismatch" if residuals["1"]["median"] > 0.05 else "radius-check-passed",
                note="Nearest cylindrical sidewall residual, excluding member ends; "
                     "not full CAD equivalence. Mesh intersections/caps can increase residuals. "
                     "Only published diameters are used by the engineering model.")
