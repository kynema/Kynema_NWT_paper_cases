"""Optional adapter for the documented GNU-Fortran Grid2Grid module interface."""

import ctypes as ct
from pathlib import Path

import numpy as np

from model import finite, wet_geometry


class Grid2Grid:
    """Interface target: LHEEA/Grid2Grid commit f5f997a35b940dd6c479b88019c2087d98d6144d."""

    def __init__(self, library, modes, zmax=16, nz=100):
        modes = Path(modes).resolve(strict=True)
        if len(str(modes).encode()) >= 300:
            raise ValueError("Grid2Grid requires a modes path shorter than 300 bytes")
        if not np.isfinite(zmax) or zmax <= 0 or nz < 4:
            raise ValueError("zmax must be positive and nz at least four")
        self.library = ct.CDLL(str(Path(library).resolve(strict=True)))
        integer, real = ct.POINTER(ct.c_int), ct.POINTER(ct.c_double)

        def bind(name, args):
            function = getattr(self.library, "__modgrid2grid_MOD_" + name)
            function.argtypes, function.restype = args, None
            return function

        initialize = bind("initializegrid2grid", [ct.c_char_p, ct.c_char_p,
                          real, real, integer, integer, real, real, integer])
        self.correct = bind("correctgrid2grid", [integer, real])
        self.get_eta = bind("gethoseta", [integer] + [real] * 4)
        self.get_u = bind("gethosu", [integer] + [real] * 7)
        end_time = bind("gethosendtime", [integer, real])
        get_depth = bind("gethoswaterdepth", [integer, real])
        self.index = ct.c_int(-1)
        values = [ct.c_double(-60), ct.c_double(zmax), ct.c_int(nz), ct.c_int(nz),
                  ct.c_double(3), ct.c_double(3)]
        initialize(b"Ocean".ljust(300, b" "), str(modes).encode().ljust(300, b" "),
                   *(ct.byref(value) for value in values), ct.byref(self.index))
        if self.index.value < 0:
            raise RuntimeError("Grid2Grid initialization did not return a valid index")
        end, depth = ct.c_double(np.nan), ct.c_double(np.nan)
        end_time(ct.byref(self.index), ct.byref(end))
        get_depth(ct.byref(self.index), ct.byref(depth))
        self.end = float(finite("HOS end time", end.value))
        self.zmax = zmax
        if self.end <= 0:
            raise ValueError("HOS record end time must be positive")
        if not np.isclose(depth.value, 60, rtol=0, atol=1e-6):
            raise ValueError(f"Expected 60 m demonstration precursor; got {depth.value} m")

    def update(self, time):
        time = float(finite("HOS time", time, ()))
        if not 0 <= time <= self.end:
            raise ValueError("Requested time lies outside the HOS record")
        value = ct.c_double(time)
        self.correct(ct.byref(self.index), ct.byref(value))

    def elevation(self, xy, time):
        xy = finite("HOS probe xy", xy)
        if xy.ndim != 2 or xy.shape[1] != 2:
            raise ValueError("HOS elevation probes must have shape (n,2)")
        self.update(time)
        result = np.empty(len(xy))
        for i, (x, y) in enumerate(xy):
            values = [ct.c_double(v) for v in (x, y, time, np.nan)]
            self.get_eta(ct.byref(self.index), *(ct.byref(v) for v in values))
            result[i] = values[-1].value
        return finite("HOS elevation", result)

    def velocity(self, xyz, time):
        xyz = finite("HOS probe xyz", xyz)
        if xyz.ndim != 2 or xyz.shape[1] != 3:
            raise ValueError("HOS velocity probes must have shape (n,3)")
        if np.any(xyz[:, 2] < -60) or np.any(xyz[:, 2] > self.zmax):
            raise ValueError("Velocity probe lies outside the Grid2Grid vertical interpolation grid")
        self.update(time)
        result = np.empty((len(xyz), 3))
        for i, (x, y, z) in enumerate(xyz):
            values = [ct.c_double(v) for v in (x, y, z, time, np.nan, np.nan, np.nan)]
            self.get_u(ct.byref(self.index), *(ct.byref(v) for v in values))
            result[i] = [v.value for v in values[-3:]]
        return finite("HOS velocity", result)


def sample(provider, geometry, time, xy_shift, time_shift=0, derivative_step=0.02):
    """Eulerian derivative at fixed wet-midpoint coordinates, not moving-node du/dt."""
    if not np.isfinite(derivative_step) or derivative_step <= 0:
        raise ValueError("derivative_step must be positive")
    time = float(finite("time", time)) + float(finite("time_shift", time_shift))
    if time - derivative_step < 0 or time + derivative_step > provider.end:
        raise ValueError("Requested time/derivative stencil lies outside the HOS record")
    shift = finite("xy_shift", xy_shift, (2,))
    endpoints = np.stack((geometry["a"], geometry["b"]), axis=1)
    eta = provider.elevation(endpoints[:, :, :2].reshape(-1, 2) + shift, time)
    eta = eta.reshape(-1, 2)
    xyz, length = wet_geometry(geometry, eta)
    wet = length > 0
    points = xyz[wet].copy()
    points[:, :2] += shift
    u, acceleration = np.zeros_like(xyz), np.zeros_like(xyz)
    u[wet] = provider.velocity(points, time)
    before = provider.velocity(points, time - derivative_step)
    after = provider.velocity(points, time + derivative_step)
    acceleration[wet] = (after - before) / (2 * derivative_step)
    return eta, u, acceleration
