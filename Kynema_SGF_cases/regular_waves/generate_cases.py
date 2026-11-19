"""Generate self-contained regular-wave mesh, aspect-ratio, and CFL sweeps."""

import argparse
import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
WAVELENGTH = 4.0
DEPTH = 4.0
LENGTH = 48.0
HEIGHT = 6.0
MAX_LEVEL = 2
BASELINE_CFL = 0.8
# Fluid properties (20 C, 1 atm); viscosities are dynamic [Pa s].
WATER_DENSITY = 1000.0
AIR_DENSITY = 1.204
WATER_VISCOSITY = 1.002e-3
AIR_VISCOSITY = 1.81e-5
# Checkpoints for restarting interrupted long runs: every CHECKPOINT_TIME_S of
# simulated time, converted to a step interval (fixed dt).
CHECKPOINT_TIME_S = 32.0
ASPECT_RATIOS = (1, 2, 4, 8)
CFLS = (0.1, 0.2, 0.5, 0.8, 1.0)
CPH_TARGETS = (1, 2.5, 5, 10, 20)
FAMILIES = {
    "linear": {
        "wave_type": "LinearWaves",
        "wave_height_m": WAVELENGTH * 0.05 / math.pi,
        "spacings": (0.0625, 0.025, 0.0125, 0.00625, 0.003125),
        # Measured/estimated CFL ratio: the orbital-velocity estimate ignores the
        # faster air motion near the interface, so dt = CFL / (factor * estimate)
        # targets the solver-reported maximum CFL. Calibrated per family.
        "cfl_rate_factor": 8.5,
        "cfl_rate_factor_source": (
            "linear/mesh_cph1, kynema-sgf bdb43d8c, realistic properties, ny=2, amg bottom "
            "solvers, fixed dt 0.0081237 s: max reported CFL 0.4026 (t=103.75 s) / estimate "
            "0.047619 = 8.455 over 0-192 s (rounded up to 8.5)"
        ),
    },
    "stokes5": {
        "wave_type": "StokesWaves",
        "wave_height_m": 0.3585,
        "spacings": (0.375, 0.15, 0.075, 0.0375, 0.01875),
        "cfl_rate_factor": 4.2,
        "cfl_rate_factor_source": (
            "stokes5/mesh_cph1, kynema-sgf bdb43d8c, realistic properties, fixed dt 0.028854 s: "
            "max reported CFL 0.7857 / estimate 0.18720 = 4.197 over 0-192 s (rounded up to 4.2)"
        ),
    },
}


def tag(value):
    return f"{value:g}".replace(".", "p")


def surface_velocities(family, count=8192):
    """Planning velocities: finite-depth Airy or deep-water Fenton fifth order.

    At kd=2*pi the deep-water approximation is adequate for timestep planning;
    accuracy analysis instead samples the solver's full finite-depth target.
    """
    a = FAMILIES[family]["wave_height_m"] / 2
    k = 2 * math.pi / WAVELENGTH
    eps = k * a
    omega = math.sqrt(9.81 * k * math.tanh(k * DEPTH))
    result = []
    for index in range(count):
        phase = 2 * math.pi * index / count
        if family == "linear":
            eta = a * math.cos(phase)
            u = a * omega * math.cosh(k * (DEPTH + eta)) / math.sinh(k * DEPTH) * math.cos(phase)
            w = a * omega * math.sinh(k * (DEPTH + eta)) / math.sinh(k * DEPTH) * math.sin(phase)
        else:
            # Deep-water limits of Fenton's H-based expansion, through order 5.
            eta = (
                eps * math.cos(phase) + eps**2 / 2 * math.cos(2 * phase)
                - 3 * eps**3 / 8 * (math.cos(phase) - math.cos(3 * phase))
                + eps**4 / 3 * (math.cos(2 * phase) + math.cos(4 * phase))
                + eps**5 * (-211 / 192 * math.cos(phase)
                            + 99 / 128 * math.cos(3 * phase)
                            + 125 / 384 * math.cos(5 * phase))
            ) / k
            amplitudes = (
                (eps - eps**3 / 2 - 37 * eps**5 / 24) * math.exp(k * eta),
                eps**4 * math.exp(2 * k * eta),
                eps**5 / 4 * math.exp(3 * k * eta),
            )
            scale = math.sqrt(9.81 / k)
            u = scale * sum(value * math.cos(n * phase) for n, value in enumerate(amplitudes, 1))
            w = scale * sum(value * math.sin(n * phase) for n, value in enumerate(amplitudes, 1))
        result.append((u, w))
    return result


def case_specs():
    cases = []
    for family, wave in FAMILIES.items():
        variants = [
            (f"mesh_cph{tag(cph)}", dx, 1, BASELINE_CFL,
             ["mesh", "aspect_ratio", "cfl"] if cph == 5 else ["mesh"], cph)
            for cph, dx in zip(CPH_TARGETS, wave["spacings"])
        ]
        baseline_dx = wave["spacings"][2]
        variants.extend(
            (f"aspect_ratio_{ar}", baseline_dx, ar, BASELINE_CFL, ["aspect_ratio"], 5)
            for ar in ASPECT_RATIOS if ar != 1
        )
        variants.extend(
            (f"cfl_{tag(cfl)}", baseline_dx, 1, cfl, ["cfl"], 5)
            for cfl in CFLS if cfl != BASELINE_CFL
        )
        k = 2 * math.pi / WAVELENGTH
        velocities = surface_velocities(family)
        for name, dx, ar, cfl, sweeps, cph in variants:
            dz = dx / ar
            coarse_dx, coarse_dz = dx * 2**MAX_LEVEL, dz * 2**MAX_LEVEL
            nx, nz = round(LENGTH / coarse_dx), round(HEIGHT / coarse_dz)
            if not math.isclose(nx * coarse_dx, LENGTH) or not math.isclose(nz * coarse_dz, HEIGHT):
                raise ValueError(f"{family}/{name}: grid does not tile the domain")
            if nx % 8 or nz % 2:
                raise ValueError(f"{family}/{name}: grid violates mesh blocking factors")
            # Kynema-SGF's VOF convention doubles the summed interface rate.
            rate = 2 * max(abs(u) / dx + abs(w) / dz for u, w in velocities)
            dt = cfl / (wave["cfl_rate_factor"] * rate)
            interval = max(1, round(0.04 / dt))
            checkpoint_interval = max(1, round(CHECKPOINT_TIME_S / dt))
            cases.append({
                "id": f"{family}/{name}",
                "family": family,
                "sweeps": sweeps,
                "wave_type": wave["wave_type"],
                "wave_height_m": wave["wave_height_m"],
                "wave_length_m": WAVELENGTH,
                "water_depth_m": DEPTH,
                "ka": k * wave["wave_height_m"] / 2,
                "domain_lo_m": [0.0, -coarse_dx, -DEPTH],
                "domain_hi_m": [LENGTH, coarse_dx, HEIGHT - DEPTH],
                "n_cell": [nx, 2, nz],
                "max_level": MAX_LEVEL,
                "dx_m": dx,
                "dy_m": dx,
                "dz_m": dz,
                "aspect_ratio": ar,
                "mesh_target_cells_per_height": cph,
                "cells_per_height": wave["wave_height_m"] / dz,
                "nominal_advective_cfl": cfl,
                "cfl_rate_estimate_per_s": rate,
                "cfl_estimate_model": "2*max_phase(|u_surface|/dx + |w_surface|/dz)",
                "cfl_rate_factor": wave["cfl_rate_factor"],
                "cfl_rate_factor_source": wave["cfl_rate_factor_source"],
                "dt_model": "nominal_advective_cfl / (cfl_rate_factor * cfl_rate_estimate_per_s)",
                "water_density_kg_m3": WATER_DENSITY,
                "air_density_kg_m3": AIR_DENSITY,
                "water_dynamic_viscosity_Pa_s": WATER_VISCOSITY,
                "air_dynamic_viscosity_Pa_s": AIR_VISCOSITY,
                "u_peak_estimate_m_per_s": max(abs(u) for u, _ in velocities),
                "w_peak_estimate_m_per_s": max(abs(w) for _, w in velocities),
                "dt_s": dt,
                "stop_time_s": 192.0,
                "output_interval_steps": interval,
                "checkpoint_interval_steps": checkpoint_interval,
                "checkpoint_interval_s": checkpoint_interval * dt,
                "sample_dt_s": interval * dt,
                "sample_x_m": [8.0, 36.0],
                "sample_points": 561,
                "analysis_time_s": [112.0, 176.0],
                "frequency_band_hz": [0.4, 0.85],
            })
    return cases


def numbers(values):
    return " ".join(f"{value:.16g}" for value in values)


def render_input(case):
    dt = f"{case['dt_s']:.16g}"
    cfl_setting = (
        f"time.cfl = {case['nominal_advective_cfl']:g}\n"
        if case["nominal_advective_cfl"] != 1.0 else ""
    )
    stokes_order = "OceanWaves.Wave1.order = 5\n" if case["family"] == "stokes5" else ""
    text = f"""# Generated by generate_cases.py; see case.json for sweep and sampling metadata.
# dt = CFL / (factor * orbital-velocity rate estimate); verify actual CFL in solver logs.
time.stop_time = {case['stop_time_s']:g}
time.fixed_dt = {dt}
{cfl_setting}time.use_force_cfl = false
time.checkpoint_interval = {case['checkpoint_interval_steps']}

transport.model = TwoPhaseTransport
transport.viscosity_fluid1 = {WATER_VISCOSITY:g}
transport.viscosity_fluid2 = {AIR_VISCOSITY:g}
transport.laminar_prandtl = 0.7
transport.turbulent_prandtl = 0.3333
turbulence.model = Laminar
incflo.physics = MultiPhase OceanWaves
incflo.gravity = 0 0 -9.81
MultiPhase.density_fluid1 = {WATER_DENSITY:g}
MultiPhase.density_fluid2 = {AIR_DENSITY:g}
ICNS.source_terms = GravityForcing
ICNS.use_perturb_pressure = true
nodal_proj.bottom_atol = 1e-14
mac_proj.bottom_atol = 1e-14
nodal_proj.bottom_solver = amg
mac_proj.bottom_solver = amg

OceanWaves.label = Wave1
OceanWaves.Wave1.type = {case['wave_type']}
{stokes_order}OceanWaves.Wave1.wave_height = {case['wave_height_m']:.16g}
OceanWaves.Wave1.wave_length = {case['wave_length_m']:g}
OceanWaves.Wave1.water_depth = {case['water_depth_m']:g}
OceanWaves.Wave1.zero_sea_level = 0.
OceanWaves.Wave1.wave_phase_offset_radians = 0.
OceanWaves.Wave1.relax_zone_gen_length = 4.
OceanWaves.Wave1.numerical_beach_length = {case['wave_length_m']:g}
OceanWaves.Wave1.numerical_beach_length_factor = 2.
OceanWaves.Wave1.initialize_wave_field = true
OceanWaves.Wave1.timeramp_period = 0.

geometry.prob_lo = {numbers(case['domain_lo_m'])}
geometry.prob_hi = {numbers(case['domain_hi_m'])}
geometry.is_periodic = 0 1 0
xlo.type = "wave_generation"
xlo.density = {AIR_DENSITY:g}
xhi.type = "pressure_outflow"
zlo.type = "slip_wall"
zhi.type = "slip_wall"

amr.n_cell = {' '.join(map(str, case['n_cell']))}
amr.max_level = {case['max_level']}
amr.blocking_factor_y = 2
amr.blocking_factor_z = 2
tagging.labels = static
tagging.static.type = CartBoxRefinement
tagging.static.static_refinement_def = static_box.txt

incflo.post_processing = sampling_fs_domain sampling_fs_target
"""
    for group, label, phase_field in (
        ("sampling_fs_domain", "fs_domain", "vof"),
        ("sampling_fs_target", "fs_target", "ow_vof"),
    ):
        frequency_setting = (
            f"{group}.output_frequency = {case['output_interval_steps']}\n"
            if case["output_interval_steps"] != 10 else ""
        )
        phase_setting = (
            f"{group}.{label}.phase_field = {phase_field}\n"
            if phase_field != "vof" else ""
        )
        text += f"""
{frequency_setting}{group}.output_format = ascii
{group}.labels = {label}
{group}.fields = {phase_field}
{group}.{label}.type = FreeSurfaceSampler
{phase_setting}{group}.{label}.plane_num_points = {case['sample_points']} 1
{group}.{label}.plane_start = 8.0 0.0 0.0
{group}.{label}.plane_end = 36.0 0.0 0.0
{group}.{label}.max_sample_points_per_cell = 32
"""
    return text + "\nincflo.verbose = 1\n"


def render_refinement(case):
    lo, hi = case["domain_lo_m"], case["domain_hi_m"]
    h = case["wave_height_m"]
    lines = [str(MAX_LEVEL)]
    for half_height in (4 * h, 2 * h):
        lines += ["1", numbers([lo[0], lo[1], -half_height, hi[0], hi[1], half_height])]
    return "\n".join(lines) + "\n"


def generate(output):
    output = Path(output)
    cases = case_specs()
    for case in cases:
        directory = output / case["id"]
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "input.inp").write_text(render_input(case))
        (directory / "static_box.txt").write_text(render_refinement(case))
        (directory / "case.json").write_text(json.dumps(case, indent=2, allow_nan=False) + "\n")
    columns = [
        "id", "cells_per_height", "dx_m", "dy_m", "dz_m", "aspect_ratio",
        "nominal_advective_cfl", "dt_s", "sample_dt_s",
    ]
    with (output / "sweep_summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(cases)
    return cases


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "input_files")
    args = parser.parse_args(argv)
    cases = generate(args.output)
    print(f"Generated {len(cases)} cases in {args.output}")


if __name__ == "__main__":
    main()
