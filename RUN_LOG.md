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
| linear/mesh_cph1 | 0.016056 | 11958 | 1993 | 0.12M | 1x36 | short | 2:00:00 | 0.4 h | 20002564 | PENDING | 20002574 | PENDING |
| linear/mesh_cph2p5 | 0.0064225 | 29895 | 4982 | 0.31M | 1x72 | standard | 6:00:00 | 1.3 h | 20002565 | PENDING | 20002575 | PENDING |
| linear/mesh_cph5 | 0.0032113 | 59790 | 9965 | 1.15M | 2x144 | hbw | 12:00:00 | 4.7 h | 20002566 | PENDING | 20002576 | PENDING |
| linear/mesh_cph10 | 0.0016056 | 119579 | 19930 | 4.3M | 8x576 | hbw | 1-00:00:00 | 8.9 h | 20002567 | PENDING | 20002577 | PENDING |
| linear/mesh_cph20 | 0.00080281 | 239159 | 39860 | 17.1M | 16x1152 | hbwl | 4-00:00:00 | 35 h | 20002568 | PENDING | 20002578 | PENDING |
| stokes5/mesh_cph1 | 0.02936 | 6540 | 1090 | 0.013M | 1x4 | short | 1:00:00 | 0.2 h | 20002569 | PENDING | 20002579 | PENDING |
| stokes5/mesh_cph2p5 | 0.011744 | 16349 | 2725 | 0.06M | 1x16 | short | 2:00:00 | 0.6 h | 20002570 | PENDING | 20002580 | PENDING |
| stokes5/mesh_cph5 | 0.0058719 | 32698 | 5450 | 0.17M | 1x36 | short | 4:00:00 | 1.5 h | 20002571 | PENDING | 20002581 | PENDING |
| stokes5/mesh_cph10 | 0.002936 | 65396 | 10899 | 0.6M | 1x72 | standard | 12:00:00 | 5.4 h | 20002572 | PENDING | 20002582 | PENDING |
| stokes5/mesh_cph20 | 0.001468 | 130792 | 21799 | 2.3M | 4x288 | hbw | 1-00:00:00 | 10.4 h | 20002573 | PENDING | 20002583 | PENDING |
