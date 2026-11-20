# Regular-wave mesh sweeps: run log

Task: run linear + stokes5 `mesh_cph*` cases (10) as set A (default) and set B
(`OceanWaves.Wave1.numerical_beach_disable_vof_term = true`), analyse, then
rebuild on `numerical_beach_upgrade_again` and rerun the better setting (set C).

## Locations (cluster kl3)

| What | Path |
|---|---|
| Case source (unmodified) | `~/Kynema_NWT_paper_cases/Kynema_SGF_cases/regular_waves/input_files/{linear,stokes5}/mesh_cph*` (repo @ 7285e77) |
| Spack env | `/scratch/mkuhn/kynema-manager/environments/env_kynema-sgf` (spec `kynema-sgf@main+hypre+mpi+netcdf+waves2amr`, develop, hash 5a77v5q) |
| Install prefix | `/kfs3/scratch/mkuhn/kynema-manager/spack/opt/spack/linux-x86_64/kynema-sgf-main-5a77v5qqe6322addq4lqncryduh52r7g` |
| Job scripts / submit dirs | `~/testruns_data/regular_waves_sweeps/` |
| Run dirs | `/scratch/mkuhn/regular_waves_sweeps/` |
| Example scripts used as template | `~/testruns_data/stokes_waves_hybrid/paper_runs/background/max_lev*/submit.sh` (account hfm, 72 ranks/node, cray FI_* env, srun block:block rank_ldom) |

Note: examples activate `env_kynema-driver`; these runs use `env_kynema-sgf`
(different develop prefix j7noqba vs 5a77v5q, so builds here do not touch
running driver-env job 19985044 `sw_sm4iw`).

## Build 1 (old code, for sets A/B)

- Source: `env_kynema-sgf/kynema-sgf`, detached HEAD at `mbkuhn/numerical_upgrade_beach_only`
- Commit: `e216ca02328cc17b7a291d9f677173c317313da0` ("enable selective use of vof term in numerical beach")
- Submodules: amrex 1fa9e823, AMReX-Hydro fa061e74, Waves2AMR 3391cd36
- Build job 19997230 (debug, 1 node) — **FAILED** (compile error):
  `relaxation_zones_ops.cpp:278: use of undeclared identifier 'disable_beach_vof_term'`.
  e216ca02 adds the local `const auto disable_beach_vof_term = wdata.disable_beach_vof_term;`
  in `modify_target_fields_for_beach()` instead of `apply_relaxation_zones()` where it is used.
  The same commit on `numerical_beach_upgrade_again` (8a39d1c9) declares it in the right function.
  Installed executable unchanged (Nov 11 build).
- User fixed and pushed: **`bdb43d8c644a627040f01e4bced11525fe1750ef`** "whoops - move declaration"
  (detached at `mbkuhn/numerical_upgrade_beach_only`, clean tree, same submodules).
- Build job 20000808 (`~/testruns_data/regular_waves_sweeps/build_bdb43d8c/build.sh`, debug 1 node):
  **COMPLETED** in 1m23s (incremental dev build reusing the spack stage dir).
- Executable: `/kfs3/scratch/mkuhn/kynema-manager/spack/opt/spack/linux-x86_64/kynema-sgf-main-5a77v5qqe6322addq4lqncryduh52r7g/bin/kynema_sgf`
  (11682184 B, 2026-11-18 12:39:17; contains `numerical_beach_disable_vof_term`).
- Archived copy: `/scratch/mkuhn/regular_waves_sweeps/exe/bdb43d8c/` (bin, lib, lib64, COMMIT).
  Binary uses DT_RPATH into the live prefix (no patchelf available), so the copy is an
  archive/identity reference, not standalone-runnable.

## Tooling

- Run-set generator: `~/testruns_data/regular_waves_sweeps/make_runs.py SET EXE_DIR [--override ...] [--tag]`
  (refuses to overwrite non-empty submit/run dirs; job refuses if run.log exists).
- Jobs `spack load kynema-sgf` from env_kynema-sgf and **abort unless the binary is
  byte-identical to the archived copy** for their build. Rebuild for step 6 only after
  all A/B jobs have left the queue.
- Analysis on cluster: `module load python/3.12.12` (numpy 2.5.2, matplotlib 3.10.8).

## Runs

Submitted 2026-11-18 (user approved plan; multi-node jobs moved to `hbw` per user).
Job IDs also in `~/testruns_data/regular_waves_sweeps/submitted_AB.txt`.
- Set A `A_default`: submit `~/testruns_data/regular_waves_sweeps/A_default/<fam>/<case>/submit.sh`,
  run `/scratch/mkuhn/regular_waves_sweeps/A_default/<fam>/<case>/` (solver log: `run.log`)
- Set B `B_novof`: same layout, deck + `OceanWaves.Wave1.numerical_beach_disable_vof_term = true`

| Case | Cells (est.) | Steps | Nodes | Ranks | Partition | Walltime | Set A job | A state | Set B job | B state |
|---|---:|---:|---:|---:|---|---|---|---|---|---|
| linear/mesh_cph1 | 0.048M | 2251 | 1 | 16 | short | 1:00:00 | 20000912 | **FAILED** (NaN, step 2) | 20000922 | HELD |
| linear/mesh_cph2p5 | 0.27M | 5627 | 1 | 36 | short | 4:00:00 | 20000913 | HELD | 20000923 | HELD |
| linear/mesh_cph5 | 1.04M | 11255 | 1 | 72 | standard | 8:00:00 | 20000914 | HELD | 20000924 | HELD |
| linear/mesh_cph10 | 4.07M | 22509 | 3 | 216 | hbw | 16:00:00 | 20000915 | HELD | 20000925 | HELD |
| linear/mesh_cph20 | 16.3M | 45018 | 10 | 720 | hbw | 2-00:00:00 | 20000916 | HELD | 20000926 | HELD |
| stokes5/mesh_cph1 | 0.005M | 2491 | 1 | 4 | short | 1:00:00 | 20000917 | **FAILED** (NaN, ~step 83) | 20000927 | HELD |
| stokes5/mesh_cph2p5 | 0.034M | 6228 | 1 | 16 | short | 1:00:00 | 20000918 | HELD | 20000928 | HELD |
| stokes5/mesh_cph5 | 0.13M | 12456 | 1 | 36 | short | 4:00:00 | 20000919 | HELD | 20000929 | HELD |
| stokes5/mesh_cph10 | 0.53M | 24913 | 1 | 72 | standard | 8:00:00 | 20000920 | HELD | 20000930 | HELD |
| stokes5/mesh_cph20 | 2.07M | 49826 | 2 | 144 | hbw | 1-00:00:00 | 20000921 | HELD | 20000931 | HELD |

Note: sacct reports these as COMPLETED because the script's final `echo` masks srun's
exit code. Check `srun exit status:` in the `%x.%j.out` file (134 = abort).

## Failures

- **A linear/mesh_cph1 (20000912)**: srun exit 134. Reported CFL 0.588 at step 1 (nominal 0.5,
  `fixed_dt does not satisfy CFL` warning), then CFL 2.99e9 at step 2, then
  `MLMG: rhs or initial residual is not finite`. Blew up immediately.
- **A stokes5/mesh_cph1 (20000917)**: srun exit 134. CFL 0.81 at steps 1-2 (nominal 0.5), then
  oscillating 1-3 (peak 24 at step 21) and back to ~0.7, then the same MLMG non-finite abort
  after ~83 steps.
- dt not changed, per README ("do not silently reduce dt in a failed run"). Asked user how to proceed.

### Solver convergence in the failed cph1 runs (user asked)

- linear/mesh_cph1: **nodal projection hit the 200-iteration cap on every solve from the
  initial projection onward** (initial proj. residual 2.32 -> 6.6e-3; step 1: 1.4e10 -> 10.2).
  MAC converged each time (24-35 iters) but its initial residual grew 7 -> 1.7e5 from the
  unconverged nodal pressure; step 2 MAC initial residual 2e71. Nodal non-convergence at t=0.
- stokes5/mesh_cph1: MAC (15-20 iters) and nodal (12-22 iters) **converged to ~1e-10 every
  step through step 81**. The nodal initial residual grew 110 -> 1617 -> 5092 -> 72233 over steps
  78-81, then step 82 MAC initial residual 1.2e14 and nodal hit the cap. Flow went unstable first;
  the solvers failed only at the end. Reported CFL was erratic from the start (0.81, then 1-3, peak 24).

### 2026-11-18: all 18 pending A/B jobs on hold (user request)

`scontrol hold` applied to 20000913-16, 20000918-31 until the failing small cases are fixed.
Release with `scontrol release <ids>`. Their submit scripts/decks are unchanged; if decks
change, cancel and restage instead of releasing.

## Test: stokes5/mesh_cph1 with dt / 10 (user request)

Reason (user): the generator's CFL estimate uses orbital velocities, but entrained air
can move faster, so dt is too large.
- Submit: `~/testruns_data/regular_waves_sweeps/tests/A_default_dtdiv10/stokes5/mesh_cph1/`
- Run: `/scratch/mkuhn/regular_waves_sweeps/tests/A_default_dtdiv10/stokes5/mesh_cph1/`
- Deck vs set A: `time.fixed_dt` 0.07706874440072348 -> 0.007706874440072347;
  sampler `output_frequency` 1 -> 10 (same 0.077 s cadence). `time.cfl = 0.5` kept (warning threshold only).
  case.json: dt_s, output_interval_steps=10, nominal_advective_cfl=0.05, test_note.
- Job 20001215 (short, 2 h): cancelled before start (user asked for debug queue).
- Job 20001243 (debug): **cancelled by user request** after 4:32 at step 4311 (t=33.2 s), stable. Max CFL 0.145 (step 3809); last-500-step mean 0.100. Partial output kept.

### dt/10 result so far (job 20001243, at ~635 steps / t~4.9 s)
Reported CFL 0.081 at step 1, max 0.122 (step 415), last 100 steps mean 0.105 / max 0.121, i.e.
~1.6-2.4x the 0.05 estimate. Projections <= 20 iters, none hit cap, no non-finite errors.

## Test: stokes5/mesh_cph1 with dt x 0.3 (= 3 x the dt/10 test; user request)
- Submit: `~/testruns_data/regular_waves_sweeps/tests/A_default_dt0p3/stokes5/mesh_cph1/`
- Run: `/scratch/mkuhn/regular_waves_sweeps/tests/A_default_dt0p3/stokes5/mesh_cph1/`
- `time.fixed_dt` 0.07706874440072348 -> 0.023120623320217044 (8304 steps); sampler
  `output_frequency` 1 -> 3 (cadence 0.0694 s vs 0.0771 s; 10/3 not an integer).
  case.json: dt_s, output_interval_steps=3, nominal_advective_cfl=0.15, test_note.
- Job 20001267 (debug): **cancelled by user request** at step 2836 (t=65.6 s), stable. Max CFL 0.6811 at step 2550 (t~59 s). 4-s binned max ~0.39 until t~24 s, then 0.51-0.68; binned mean rose 0.27 -> ~0.36 (not plateaued). Partial output kept.

## Test: stokes5/mesh_cph1 with dt targeting measured CFL 0.85 (user request)
- Rate = highest measured CFL per unit dt = 0.6811 / 0.0231206 = 29.46 /s (dt x0.3 test, t~59 s).
  (dt/10 test gave 18.77 /s but only reached t=33 s.)
- `time.fixed_dt` = 0.85 / 29.46 = **0.02885410339477975 s** (0.374 x generated dt; 6654 steps;
  orbital-estimate nominal CFL 0.187). Sampler every step (0.0289 s).
- Submit: `~/testruns_data/regular_waves_sweeps/tests/A_default_cfl0p85/stokes5/mesh_cph1/`
- Run: `/scratch/mkuhn/regular_waves_sweeps/tests/A_default_cfl0p85/stokes5/mesh_cph1/`
- Job 20001301 (debug): **COMPLETED**, srun exit 0, 8:32 wall, reached t=192.02 s (6655 steps).
  CFL: step 1 0.304; max **0.782** (step 4800, t~138 s); never > 0.85. 16-s binned max 0.50 (0-16 s),
  0.64, then 0.67-0.78 from t=32 s on; binned mean 0.34 -> ~0.42-0.43 from t~64 s (plateaued).
  Projections <= 20 iters, no cap hits, no non-finite. 26624 sampler files (both groups, every step).

## Realistic properties + CFL factor (user request, 2026-11-18)

User: air undergoes strong rotational motion, so use realistic viscosities and air density;
calibrate a CFL factor from the debug case; all non-CFL-sweep cases target CFL 0.8 (baseline
shared with CFL sweep, user chose 0.8 not 0.85); ignore diffusion in dt; drop `incflo.diffusion_type`.

- Test `A_default_visc_cfl0p85` (job 20001385): **FAILED** at startup, density check
  (rho_g_max 1.204 vs rho_g_min 1.0).
- Test `A_default_visc2_cfl0p85` (job 20001736): same failure with zlo/zhi zero_gradient density/vof BCs.
- Cause: `OceanWavesBoundary::set_density` blends rho_water*target_vof with the existing Dirichlet
  ghost value, assumed to be gas density = `xlo.density` (default 1.0). Fix: `xlo.density = 1.204`.
- Test `A_default_visc3_cfl0p85` (job 20001761, debug): **COMPLETED**, srun exit 0, 12:22 wall, t=192 s.
  Properties: water rho 1000, mu 1.002e-3; air rho 1.204, mu 1.81e-5; default diffusion; xlo.density 1.204.
  dt 0.02885410339477975 (same as inviscid test). Max CFL **0.7857** (t=186 s), conv only (diff 0, implicit);
  16-s binned max 0.49 -> 0.62-0.79, mean 0.33 -> ~0.41. Projections <= 20 iters, no cap, no non-finite.
  (Inviscid run at same dt: max 0.7816.) Actual grid 12544 cells; ~28k cell-steps/s/core.
- **CFL factor F = 0.7857 / (dt * 6.4877/s = 0.18720) = 4.197 -> 4.2.** dt = C / (F * orbital estimate).

### Generator/repo changes (local working tree, uncommitted)
`generate_cases.py`: BASELINE_CFL 0.8; properties; `xlo.density`; no diffusion_type; CFL_RATE_FACTOR 4.2
and metadata (cfl_rate_factor, source, dt_model, properties). CFL sweep now cfl_0p1/0p2/0p5/1
(+ mesh_cph5 at 0.8); removed stale `input_files/*/cfl_0p8`. Tests updated (21 pass, 1 plot skip).
README + input_files/README updated. All dt are x0.381 of the old values for mesh cases.

## Linear cph1 investigation (v2 decks: F=4.2, CFL 0.8, realistic properties)

- `tests/A_v2_debug` (job 20002019): nodal projection hit the 200-iter cap on all 110 solves (from t=0);
  CFL 0.22 -> 1.5 at step 4, 1.5-3.9 after; NaN at step 198 (t=6.4 s), srun 134.
- `tests/A_v2_semicoarsen_debug` (20002142, + nodal_proj.do_semicoarsening=1, max_semicoarsening_level=30):
  no change (all 109 nodal solves capped). Cancelled at step 136.
- `tests/A_v2_nz32_debug` (20002190, nz 24->32, top z=4): all 37 nodal solves capped. Cancelled (user meant y).
- `tests/A_v2_ny8_debug` (20002210, **ny 2->8 at same dy, span 2 m**, static boxes widened):
  **nodal converges in 8-11 iters, MAC 13-16, no cap hits.** CFL path unchanged (0.23, 0.23, 0.36, 1.53, ...),
  max 5.1; NaN at step 273 (t=8.87 s).
- `tests/A_v2_ny8_plt_debug` (20002306, 10 steps, plotfile every step; locator
  `~/testruns_data/regular_waves_sweeps/plt_cfl_locate.py`): CFL peak always in pure-air cells of the first
  level-2 layer above SWL (z=0.031 m), |u|,|w| ~0.7-1.3 m/s (orbital 0.13 m/s), at x~3.2-3.4 (gen-zone edge),
  7.2, 39-40; uniform in y. Plotfile CFL reproduces solver log.
- Stokes5 cph1 also has ~1 finest cell per H (dz 0.375 vs H 0.3585); difference vs linear not yet explained.
- User: calibrate the linear dt separately. `tests/A_v2_ny8_dtq_debug` (job **20002346**, debug, 1x72,
  55 min): ny8 deck, dt = 0.008123711632166726 (F4.2 dt / 4), sampler every 5 steps.
  Cancelled by user request at step 6812 (t=55.3 s), stable. CFL 0.13-0.14 to t~20 s, then rose;
  **max 0.4377 at t~51 s**. Nodal 8-11 iters (mean 10.0), MAC 13-20 (mean 14.7), no cap hits.

### Solver tests at ny=2 (user request): native AMG bottom solver x nsolve x semicoarsening
Linear cph1 v2 deck (ny=2), dt = 0.008123711632166726, `time.max_step = 20`, same settings on
nodal_proj and mac_proj: `bottom_solver = amg`, `do_nsolve` on/off, `do_semicoarsening = 1` +
`max_semicoarsening_level = 30` on/off. Watcher scancels on any 200-iteration solve.
Dirs `tests/solver_amg_ns{0,1}_sc{0,1}/linear/mesh_cph1`; jobs 20002391 (ns0 sc0), 20002393 (ns0 sc1),
20002394 (ns1 sc0), 20002395 (ns1 sc1).
Results (all srun exit 0, 20 steps, no cap hits; CFL identical in all four and to the ny8 dt/4 run, 0.057 -> 0.124):
| test | nodal iters mean (range) | MAC iters mean (range) | wall/step |
|---|---|---|---|
| amg | 10.3 (8-11) | 13.9 (11-15) | 0.107 s |
| amg + semicoarsening | 93.1 (63-124) | 23.2 (21-27) | 0.334 s |
| amg + nsolve | 10.3 (8-11) | 13.9 (11-15) | 0.107 s |
| amg + nsolve + semicoarsening | 93.1 (63-124) | 23.2 (21-27) | 0.334 s |
-> AMG bottom solver fixes the ny=2 nodal convergence; semicoarsening hurts; nsolve no effect.

### v3 decks: AMG bottom solver (user approved)
Generator adds `nodal_proj.bottom_solver = amg` and `mac_proj.bottom_solver = amg` to all cases;
ny stays 2. Tests (21 pass) and README updated. Local working tree, uncommitted.

### Linear calibration on v3 deck (ny=2, amg)
`tests/A_v3_amg_dtq_debug/linear/mesh_cph1`, job **20002421** (debug, 1x36, 55 min):
dt = 0.008123711632166726 (generated dt / 4), 23635 steps, sampler every 5 steps.
**Result**: COMPLETED, srun exit 0, 34:49 wall, t=192.004 s (23635 steps), no non-finite.
Nodal mean 10.0 (max 11), MAC mean 13.5 (max 16), no cap hits. CFL ~0.13-0.14 to t~40 s, rises
45-55 s, peak **0.4026 at t=103.75 s**, then drifts to ~0.36; mean ~0.26 late; max in 112-176 s 0.3885.
**Linear factor = 0.4026 / (0.0081237 * 5.8617 = 0.047619) = 8.455 -> 8.5 proposed** (stokes5 4.2).

## v3 decks: per-family CFL factor + checkpoints (user approved)
- Generator: linear F=8.5 (cph1 calibration above), stokes5 F=4.2; `time.checkpoint_interval` =
  round(32 s / dt) per case (restart via `io.restart_file`; step count continues). Tests 21 pass; README updated.
- Uploaded regenerated decks to `~/testruns_data/regular_waves_sweeps/input_files_v3/` (cluster repo copy untouched).
- **Cancelled the 18 held A/B jobs** (20000913-16, 20000918-31): old dt/properties/solver settings, never ran.
- **Submitted 2026-11-18** (user approved; IDs also in `~/testruns_data/regular_waves_sweeps/submitted_AB_v3.txt`): `A_default_v3`, `B_novof_v3` (B adds numerical_beach_disable_vof_term = true).
  Run dirs `/scratch/mkuhn/regular_waves_sweeps/{A_default_v3,B_novof_v3}/<fam>/<case>/`.

| Case | dt [s] | Steps | Ckpt every | Cells (est.) | Nodes x ranks | Partition | Walltime | Est. run | A job | A state | B job | B state |
|---|---:|---:|---:|---:|---|---|---|---:|---|---|---|---|
| linear/mesh_cph1 | 0.016056 | 11958 | 1993 | 0.12M | 1x36 | short | 2:00:00 | 0.4 h | 20002564 | COMPLETED (0) | 20002574 | COMPLETED (0) |
| linear/mesh_cph2p5 | 0.0064225 | 29895 | 4982 | 0.31M | 1x72 | standard | 6:00:00 | 1.3 h | 20002565 | COMPLETED (0) | 20002575 | COMPLETED (0) |
| linear/mesh_cph5 | 0.0032113 | 59790 | 9965 | 1.15M | 2x144 | hbw | 12:00:00 | 4.7 h | 20002566 | COMPLETED (0) | 20002576 | COMPLETED (0) |
| linear/mesh_cph10 | 0.0016056 | 119579 | 19930 | 4.3M | 8x576 | hbw | 1-00:00:00 | 8.9 h | 20002567 | CANCELLED (user; stuck sampler pts) | 20002577 | CANCELLED (user; stuck sampler pts) |
| linear/mesh_cph20 | 0.00080281 | 239159 | 39860 | 17.1M | 16x1152 | hbwl | 4-00:00:00 | 35 h | 20002568 | **NODE_FAIL** step 46021 | 20002578 | CANCELLED (user; stuck sampler pts) |
| stokes5/mesh_cph1 | 0.02936 | 6540 | 1090 | 0.013M | 1x4 | short | 1:00:00 | 0.2 h | 20002569 | COMPLETED (0) | 20002579 | COMPLETED (0) |
| stokes5/mesh_cph2p5 | 0.011744 | 16349 | 2725 | 0.06M | 1x16 | short | 2:00:00 | 0.6 h | 20002570 | COMPLETED (0) | 20002580 | COMPLETED (0) |
| stokes5/mesh_cph5 | 0.0058719 | 32698 | 5450 | 0.17M | 1x36 | short | 4:00:00 | 1.5 h | 20002571 | COMPLETED (0) | 20002581 | COMPLETED (0) |
| stokes5/mesh_cph10 | 0.002936 | 65396 | 10899 | 0.6M | 1x72 | standard | 12:00:00 | 5.4 h | 20002572 | COMPLETED (0) | 20002582 | COMPLETED (0) |
| stokes5/mesh_cph20 | 0.001468 | 130792 | 21799 | 2.3M | 4x288 | hbw | 1-00:00:00 | 10.4 h | 20002573 | CANCELLED (user; stuck sampler pts) | 20002583 | CANCELLED (user; stuck sampler pts) |

### Status 2026-11-19 (after reconnect)
14/20 COMPLETED with srun exit 0; no projection cap hits, no non-finite values in any job.
- **A linear/mesh_cph20 (20002568): NODE_FAIL** (hardware, node x1005c2s4b1n1) after 9:47:22 at step 46021
  (t=36.9 s). Checkpoints chk00000, chk39860 (3.2 GB each) intact. Observed rate ~4700 steps/h
  (full run ~51 h vs 35 h estimate).
  **Restart job 20007733** (hbwl, 16x1152, 4 days) from `.../A_default_v3/linear/mesh_cph20/chk39860`
  in a NEW run dir `/scratch/mkuhn/regular_waves_sweeps/A_default_v3/linear/mesh_cph20_r1/`
  (submit dir `~/testruns_data/.../A_default_v3/linear/mesh_cph20_r1/`), so the failed run's files are
  not overwritten. Its post_processing starts after step 39860; the 112-176 s window is fully covered.
- Max reported CFL after step 50 (complete or so far), target 0.8:
  linear A/B: cph1 0.830/0.836, cph2p5 0.774/0.704, cph5 0.350/0.345, cph10 0.196/0.190 (running),
  cph20 0.146/0.126 (partial); stokes5 A/B: cph1 0.794/0.764, cph2p5 0.351/0.350, cph5 0.297/0.297,
  cph10 0.322/0.322, cph20 0.317/0.314 (running).
  -> factors calibrated at cph1 hold at cph1 (linear also cph2p5) but over-correct at finer meshes, so the
  finer runs sit at CFL ~0.15-0.35 (stable; mesh sweep is not at constant achieved CFL).

### A/B v3 analysis attempt and FreeSurfaceSampler bug (2026-11-19)
- `compare_sets.py` on the 7 completed pairs (`/scratch/mkuhn/regular_waves_sweeps/comparison_AB_v3_partial/`):
  target apparent Cr up to 0.51 -> invalid for 4 of 7 pairs.
- Cause: **stuck FreeSurfaceSampler points** (constant -4 m = domain bottom, or the refined-band bottom edge),
  identical in A/B and in surface/target. Stuck counts (of 561): linear cph1 0, cph2p5 267, cph5 90, cph10 23,
  cph20 166; stokes5 cph1 0, cph2p5 28, cph5 25, cph10 0, cph20 2.
- **Root cause (confirmed):** `FreeSurfaceSampler::initialize` assigns points to cells with
  `ceil((xm -/+ 0.5*dx - s)/dxs)`; neighbours compute a shared face from different centres, rounding can leave a
  point on a face owned by no cell; `update_sampling_locations` initialises heights to plo, so it stays at the
  bottom (or gets the top of an uncovered coarse water cell). An icpx (-O3 -fno-math-errno) replica
  (`~/testruns_data/regular_waves_sweeps/assign_test.cpp`) reproduces all 8 nonzero counts and the exact stuck x
  positions; gcc gives different counts. Proposed fix: compute bounds from integer face indices
  (`ceil((plo + f*dx - s)/dxs)` for f=i, i+1), same for the y test; add an ownership check.
- User: **cancelled** the running high-res jobs (A/B linear cph10, B linear cph20, A/B stokes5 cph20) and the
  pending/started A linear cph20 restart 20007733 (ran 13:33).
- Clean pairs (A vs B, 112-176 s): linear cph1 Cr 0.026/0.165, L2/H 6.0%/12.4%; stokes5 cph1 Cr 0.170/0.218;
  stokes5 cph10 Cr 0.0225/0.0154, L2/H 20.3%/19.6% (within noise).
- Window scan (32 s windows): Cr steady from 32-48 s in resolved cases; stokes5 cph1 settles ~112 s; stokes5
  cph10 phase error and L2 drift monotonically (-19 -> -31 deg, 14 -> 21%) in both sets.

### v4 decks and analysis (user request)
- Sampling over the whole run, native format, ~0.08 s cadence (`round(0.08/dt)` steps); default analysis window
  **64-190 s**, halves 64-127 / 127-190. Native reader in `reflection_coeff_komega.py` (ASCII still supported);
  `analyze_case.py` rejects stuck points. Tests 24 pass. Verified on real native output
  (`tests/native_fmt_debug`, job 20008027). Code copy for cluster analysis:
  `~/testruns_data/regular_waves_sweeps/analysis_code_v4/`. Generator/README regenerated locally (uncommitted).

### kynema-sgf sampler fix: branches (user request, 2026-11-19)
- Fix verified before commit: icpx replica of the new rule gives every point exactly one owner in all 10 meshes
  (x and y), patched file compiles with the build flags with no warnings. Patch:
  `~/testruns_data/regular_waves_sweeps/sampler_fix/freesurfacesampler_owner_fix.patch`.
- New branch **`fix_free_surface_sampler_ownership` = `3ece7a3f`** from origin/main `08e13265`
  (FreeSurfaceSampler.cpp only: per-point host-computed owner index; count/store loops in initialize() and
  post_regrid_actions(); update_sampling_locations() unchanged).
- Rebased onto it (no conflicts):
  - `numerical_upgrade_beach_only`: bdb43d8c -> **`13de3174`** (content differs from bdb43d8c only by the fix).
  - `numerical_beach_upgrade_again`: e83aee23 (fork tip, newer than 8a39d1c9 seen earlier) -> **`e3f213f7`**;
    also picks up 4 origin/main commits it was behind (08e13265, cc11c499, b9fa2163, 4a1df971).
- Pushed to fork `mbkuhn` (from the local Mac via SSH; kl3 has no GitHub credentials):
  new branch + force-with-lease on both rebased branches (leases pinned to bdb43d8c / e83aee23).
- env_kynema-sgf source now on `numerical_upgrade_beach_only` @ 13de3174, submodules updated, tracking mbkuhn.
  **Not yet rebuilt.**

## Build 2 for A/B: 13de3174 (sampler fix)
- Build job 20008114 (`~/testruns_data/regular_waves_sweeps/build_13de3174/`, debug): COMPLETED 0:55, no errors.
  Executable (same spack prefix 5a77v5q) 11670008 B, 2026-11-19 10:14:17; differs from bdb43d8c archive.
  Archive: `/scratch/mkuhn/regular_waves_sweeps/exe/13de3174/` (COMMIT file); old `exe/bdb43d8c` kept.
- Sampler check `tests/v4_fix_debug/linear/mesh_cph2p5` (job 20008120, v4 deck, 300 steps): srun exit 0;
  **0 stuck points** in both sampler groups (previously 267), target z exactly +/-0.0318 m, no -4 m values.

## v4 run sets (staged, not submitted)
`A_default_v4`, `B_novof_v4` from `input_files_v4` (native output ~0.08 s, window 64-190 s, checkpoints ~32 s),
build 13de3174, same resources as v3. Generator: `make_runs.py` (SRC input_files_v4, job names rwA4_/rwB4_).

## v5: recalibrated dt (user request): per-family well-resolved F + stability multipliers
- Ratios measured/estimated max CFL from v3 (max of A/B): linear cph1 8.88, 2.5 8.22, 5 3.72, 10 2.08,
  20 >=1.55 (ended at 37 s, pre-plateau); stokes5 cph1 4.17, 2.5 1.84, 5 1.56, 10 1.69, 20 1.67.
- User chose per-family F (linear 2.1, stokes5 1.9) + multipliers (linear cph1 x4.3, cph2.5 x4.0, cph5 x1.8;
  stokes5 cph1 x2.25); linear cph20 trusts the trend, check CFL after ~55 s.
- dt (s) / steps: linear cph1 0.015114/12.7k, 2.5 0.006499/29.5k, 5 0.007221/26.6k, 10 0.006499/29.5k,
  20 0.0032495/59.1k; stokes5 cph1 0.028844/6.7k, 2.5 0.02596/7.4k, 5 0.01298/14.8k, 10 0.00649/29.6k,
  20 0.003245/59.2k. Tests 24 pass; README updated.
- v3 throughput (steps/h): lin cph1 36.3k (1n), 2.5 19.3k (1n), 5 12.2k (2n), 10 7.9k (8n), 20 4.76k (16n);
  stokes5 cph1 56.3k, 2.5 43.3k, 5 28.8k, 10 16.1k (1n), 20 11.8k (4n). Actual cells: lin 0.12/0.40/1.38/4.76/17.5M;
  stokes5 0.013/0.057/0.18/0.61/2.25M.
- Decks `~/testruns_data/regular_waves_sweeps/input_files_v5/`; staged `A_default_v5`, `B_novof_v5`
  (build 13de3174, native ~0.08 s, window 64-190 s). Linear cph20 back on hbw (16 nodes, 36 h).
  v4 staging dirs are superseded (never submitted).

### v5 submitted 2026-11-19 (user approved). IDs: `~/testruns_data/regular_waves_sweeps/submitted_AB_v5.txt`
Status snapshot: `~/testruns_data/regular_waves_sweeps/run_status.sh submitted_AB_v5.txt`

| Case | A job | A state | B job | B state |
|---|---|---|---|---|
| linear/mesh_cph1 | 20008260 | COMPLETED (0), maxCFL 0.764 | 20008271 | COMPLETED (0), maxCFL 0.740 |
| linear/mesh_cph2p5 | 20008262 | RUNNING (CFL ~0.72) | 20008272 | RUNNING (CFL ~0.73) |
| linear/mesh_cph5 | 20008263 | CANCELLED (unstable) | 20008273 | CANCELLED (unstable) |
| linear/mesh_cph10 | 20008264 | **FAILED** NaN step 2 | 20008274 | **FAILED** NaN step 2 |
| linear/mesh_cph20 | 20008265 | CANCELLED (unstable) | 20008275 | CANCELLED (unstable) |
| stokes5/mesh_cph1 | 20008266 | COMPLETED (0), maxCFL 0.776 | 20008276 | COMPLETED (0), maxCFL 0.812 |
| stokes5/mesh_cph2p5 | 20008267 | COMPLETED, **maxCFL 1.80 (unusable)** | 20008277 | COMPLETED, **maxCFL 1.80 (unusable)** |
| stokes5/mesh_cph5 | 20008268 | CANCELLED (unstable) | 20008278 | CANCELLED (unstable) |
| stokes5/mesh_cph10 | 20008269 | CANCELLED (unstable) | 20008279 | CANCELLED (unstable) |
| stokes5/mesh_cph20 | 20008270 | CANCELLED (unstable) | 20008280 | CANCELLED (unstable) |

### v5 outcome: recalibration unstable on finer meshes
- Stable near 0.8: linear cph1, stokes5 cph1 (A/B), linear cph2.5 (running ~0.72-0.73).
- Unstable (first CFL>1 at step 3 for linear cph5/10/20; t~1.3-1.4 s for stokes5 cph5/10/20; step 8 for
  stokes5 cph2.5): linear cph10 NaN at step 2; linear cph20 CFL 2.2-4.0 from start; linear cph5 steady
  1.3-1.55; stokes5 cph5/10/20 periodic bursts to 7-10. Cancelled the 10 still-running unstable jobs (watcher
  had silently failed after its first SSH error). Ruled out: gravity number dt*sqrt(g/dz), explicit diffusion.
- v3 dt (lin F 8.5, stokes 4.2) was stable for all these meshes (CFL 0.15-0.35).

### dt stability scan (user chose option 2)
`stage_dt_scan.py`: dt = r x v3 dt, set-A v5 deck, stop 10 s, partition short; `tests/dtscan/<tag>`;
IDs `~/testruns_data/regular_waves_sweeps/submitted_dtscan.txt`; status `dtscan_status.sh`.
r ladders: stokes5 cph2.5/5/10/20 and linear cph5: 1.3, 1.6, 1.9; linear cph10/20: 1.5, 2, 3.
(v5 dt corresponds to r = 2.21 stokes5, 2.25 linear cph5, 4.05 linear cph10/20.)

### dt-scan and v5 runs cancelled (user switched to adaptive dt)
Partial scan results before cancel (t reached, max CFL, first CFL>1): stokes5 cph2.5 r1.3 done (0.445, never);
stokes5 cph10 r1.3 (6.5 s, 0.425, never), r1.6 (7.9 s, 1.65, **7.6 s**), r1.9 (9.0 s, 1.80, **6.0 s**);
stokes5 cph20 r1.3/1.6 (1.6/1.9 s, 0.42/0.74, never); linear cph5 r1.3/1.6/1.9 (5.5-6.6 s, 0.32/0.40/0.76, never);
linear cph10 r1.5/2 (1.1/1.7 s, 0.24/0.44, never). All v5 linear cph2.5 A/B cancelled at ~t=95 s (CFL <=0.76).

## v6: adaptive dt (user request)
- `time.fixed_dt = -1`, `time.initial_dt` = calibrated estimate (per-family F + multipliers, = v5 dt),
  `time.cfl` = target (0.8; sweep values), always explicit. Solver: initial_dt only used if CFL~0 at t=0;
  first dt shrunk x0.1 (`init_shrink`), growth <=10%/step (`max_dt_growth`).
- Output by time: samplers `output_time_interval = 0.08`, checkpoints `time.checkpoint_time_interval = 32`.
  case.json: time_stepping=adaptive, initial_dt_s, sample_time_interval_s, checkpoint_time_interval_s
  (fixed-dt fields dropped).
- Analysis: output times read from native Header (full precision); gaps >2.5x median rejected (incl. outputs
  just outside the window); accuracy norms time-weighted; harmonic fits on actual times; k-omega on data
  linearly resampled to a uniform grid (flag in results). New non-uniform test; 25 tests pass. README updated.
- Cluster: decks `input_files_v6/`, analysis code `analysis_code_v6/`, sets `A_default_v6`, `B_novof_v6`
  (build 13de3174); walltimes cover v3-level dt (linear cph20 on hbwl, 16 nodes, 3 days).
- **Submitted** 20 jobs 20008806-20008825 (`submitted_AB_v6.txt`); status `run_status_adaptive.sh`.

### v6 progress
- Completed cleanly (CFL capped at 0.8, no solver cap hits, exit 0): linear cph1, cph2.5; stokes5 cph1, cph2.5,
  cph5, cph10 (A and B); A linear cph5. Interim A/B comparison of the 6 complete pairs:
  `/scratch/mkuhn/regular_waves_sweeps/comparison_AB_v6_partial/`.

### stokes5 adaptive-dt dips: located at the wave-generation inlet
- Observation: in stokes5 runs (cph5/10/20, both A and B at similar rates; none in linear) dt intermittently drops
  2-6x for one step, then regrows at 10%/step (~15 steps). Spike is purely convective; velocity-solve initial
  residual is 3-4x normal on the step before each dip. Dips recur at a fixed wave phase (A cph10: 32.186 s and
  36.801 s are 3 periods apart, same minimum dt 0.0014079 s).
- Plot run: `tests/A_v6_s10_plt_dtdip` (job 20010692, short, 72 ranks): restart of A_default_v6 stokes5 cph10 from
  chk04186, plotfile every step to step 4850 (density velocity vof p). 7 dips (onset steps 4220, 4396, 4409, 4450,
  4646, 4790, 4814); not bit-identical to the original after restart but same pattern (4814 at t=36.799).
  Scripts: `plt_rate_scan.py` (vectorized max-rate locator: scan/top/diff), `plt_column.py`.
- Result: in all 7 dips the fastest cell is on level 2 in the first or second cell from xlo (x = 0.019 or 0.056 m),
  in air (vof 0) just above the interface, uniform across y. Rates 76-243 /s vs ~50 /s normal (normal max is air
  just above a crest anywhere in the tank or beach). Largest: u = -6.8 m/s, w = -2.3 m/s (cell-centred; air moving toward the inlet face within the
  first cell, not flux through it: the face normal velocity is Dirichlet) at x=0.019, z=0.219, at the step where the rising surface first enters the cell below
  (vof 0 -> 0.006 at z=0.181); one step later u = -2.2, then normal. Neighbouring cells show a short-lived
  recirculation (w up to +2.9 m/s at x=0.056).
- Interpretation (unverified; ghost cells not in plotfiles): an artefact of the wave_generation boundary/relaxation
  at the air-water interface when the surface crosses a cell face at the inlet, not physics in the test section.
  Affects dt only (x<0.06 m, far upstream of the 8-36 m analysis span).
- Code check (13de3174, OceanWavesBoundary.cpp): wave_generation ghost/face velocity = full target liquid velocity
  (ow_velocity at the ghost cell centre) wherever target ow_vof > VOF_TOL (1e-12), else left at the Dirichlet
  value 0 (air). So as the target surface rises past a z-face at the inlet, the boundary face velocity of that
  row switches from 0 to the full wave velocity in one step (and back on the falling side), with no vof weighting.
  That step change in boundary flux into a light (rho 1.2) air cell is the likely trigger of the one-step
  cell-centred spike; consistent with timing (surface entering a new cell row at the inlet). Not judged a clear
  bug (design choice). Set C (e3f213f7) does not touch the boundary fill but changes ow_velocity in partially
  filled cells to the liquid-centroid height (update_target_velocity / approx_liquid_centroid_z_interp), so the
  value switched on drops from the cell-centre value (above the surface, possibly extrapolated) to between the
  cell below and the centre; it reduces but does not remove the on/off jump. No changes until C runs.
- Decision deferred to after set C (user): whether to blend air/water velocity in two-phase inlet ghost cells.
  Pro: smoother boundary velocity, consistent with how the relaxation zone blends phase velocities (which
  dominates the inlet region anyway). Con: boundary value deviates further from wave theory. Test: compare
  stokes5 dt-dip frequency/depth in C vs the matching v6 set (same counting as above: steps with dt < 0.3/0.2/0.1
  x median, and onset locations via plt_rate_scan.py if needed).

### compare_sets.py: load-once + corrected noise windows
- Found: compare_sets.py (cluster) still used early/late = 112-144 / 144-176 s for the "±" noise (the main
  metrics already used the full 64-190 s window from case.json). The ± values in the interim v6 A/B tables
  (comparison_AB_v6_partial, _partial2, _partial3) are therefore half-spreads of those two 32 s subwindows,
  not of 64-127 / 127-190 s as previously stated.
- New compare_sets.py (old kept as compare_sets_v1.py): windows results (64-190), half1 (64-127),
  half2 (127-190), written as analysis/half1|half2.{json,npz} (old early/late files untouched). Sampler file
  readers in reflection_coeff_komega (find_profiles, read_native_time, read_native_profile, read_ascii_profile)
  are memoized per case, so each case is read from disk once for all windows; analysis code unchanged.
- Check (A stokes5 cph5): cold full window 123 s, then half1 0.3 s and full again 0.7 s; recomputed summary
  matches stored results.json to ~1e-15 relative (last-digit floating-point differences from the fits, e.g.
  rel. amp. error -0.17369413524636568 vs ...658); Cr identical.

## Set C (user: code C with option A) and conditional set D
- User decision (Nov 19, 23:xx MST): run C = numerical_beach_upgrade_again with A's deck (default beach VOF
  term). If the stokes5 dt drops improve in C, run D = same code with option B
  (`numerical_beach_disable_vof_term = true`). Check dt drops early on the coarser stokes5 cases only.
- Interim nine-pair A/B result at the time (linear cph20 still running): linear favours B from cph2.5 up;
  stokes5 mixed (A cph2.5/5, B cph10, ~tie cph20).
- Build deviation from the original step 6 ("switch the env's kynema-sgf"): the two v6 linear cph20 jobs still
  run on the env_kynema-sgf install (shared libs libamrex_3d.so, libWaves2AMR.so in prefix 5a77v5q), so that
  env is NOT touched. New env `/scratch/mkuhn/kynema-manager/environments/env_kynema-sgf-C` (copied spack.yaml +
  include.yaml) with its own full clone of github.com/mbkuhn/kynema-sgf, branch numerical_beach_upgrade_again
  @ **e3f213f72784d51a39777dd7c1dcd0b00c6b58b6**; submodules amrex 1fa9e823, AMReX-Hydro fa061e74,
  Waves2AMR **ad6861d1** (A/B: 3391cd36), googletest 52204f78, doxygen-awesome-css af1d9030.
  (Cloning from the env repo failed: it is a blob:none partial clone.) Concretized: only kynema-sgf to build,
  new prefix kynema-sgf-main-**jjkukzg** (dev_path differs); hypre xnjh7cs, netcdf-c go67emj, fftw jmkj35f reused.
  All OceanWaves.Wave1.* deck keys used by A/B exist in C (incl. numerical_beach_disable_vof_term).
- Build job 20014160 (`~/testruns_data/regular_waves_sweeps/build_C_e3f213f7/build.sh`, debug, `spack install -j 104`).
- make_runs.py gained `--env` (old kept as make_runs_v6AB.py); submit scripts record and activate the env.
- Build 20014160 COMPLETED 7:04, no errors: prefix kynema-sgf-main-jjkukzg, exe 11724904 B (2026-11-19 23:10),
  links its own libamrex_3d.so / libWaves2AMR.so; contains update_target_velocity. A/B exe re-checked: still
  byte-identical to exe/13de3174. Archive: `/scratch/mkuhn/regular_waves_sweeps/exe/e3f213f7/` (bin, lib, lib64, COMMIT).
- Set C staged: `make_runs.py C_default_v6 exe/e3f213f7 --env env_kynema-sgf-C --tag e3f213f7`; decks identical
  to A_default_v6 (input_files_v6, no overrides), same resources.
- **Submitted** 10 jobs 20014199-20014208 (`submitted_C_v6.txt`, stokes5 first).
- Early dt-dip check: `dtdip_compare.sh CASE SET...` (common simulated span; steps with dt < 0.3/0.2/0.1 x median,
  dip events = dt/prev < 0.7, min dt). Reference over the full run, stokes5 cph10: A 213/48/0, 230 events;
  B 198/37/0, 224 events. Plan: judge on stokes5 cph10 at t >= 40 s (A had 9 dips < 0.3 x median by 40 s);
  cph5 too sparse early (A first dip 89 s).

### Nov 20 morning status (connection to kl3 lost overnight ~23:30-09:00 MST; jobs unaffected)
- A/B v6 linear cph20 COMPLETED cleanly (20008810, 20008820): **all 20 A/B v6 runs complete**, all clean.
- Set C: 9/10 COMPLETED cleanly (CFL capped 0.8, no caps, no nonfinite, exit 0); linear cph20 (20014208) running, t=88 s at 09:08.
- Corrected-halves nine-pair A/B comparison finished (comparison_AB_v6_partial3_halves, no errors).
- **dt dips, full runs (dtdip_compare.sh), stokes5 steps <0.3/<0.2/<0.1 x median; events:**
  cph1 A 1/0/0, B 5/0/0, C 12/2/0; cph2.5 all 0; cph5 A 29/10/0 (108), B 42/16/0 (109), C 19/6/0 (89);
  cph10 A 213/48/0 (230), B 198/37/0 (224), C 185/34/0 (208); cph20 A 853/295/9 (678), B 740/235/2 (696),
  C 851/295/10 (726). A vs B (identical inlet code) spread ~10-25% = run-to-run variability.
  C within that spread at cph10/20, better at cph5, worse (shallow) at cph1 -> **no clear improvement; D not launched**.
- Final comparisons started under nohup (`run_final_comparisons.sh`): A vs B all 10 pairs -> comparison_AB_v6_final;
  A vs C nine pairs -> comparison_AC_v6_partial; status lines in final_comparisons.status.
