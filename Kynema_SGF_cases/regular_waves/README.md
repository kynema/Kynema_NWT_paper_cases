# Regular-wave accuracy and reflection sweeps

Only **linear** and **fifth-order Stokes** waves are used. Every new run supplies
both accuracy and k-omega reflection measurements; there are no separate
reflection-only or convergence-only inputs.

## Physical setup and cost

| Parameter | Linear | Fifth-order Stokes |
|---|---:|---:|
| Wavelength, lambda | 4 m | 4 m |
| Wave height, H (crest to trough) | 0.0636619772 m | 0.3585 m |
| Steepness, ka = pi H / lambda | 0.05 | 0.281565 |
| Water depth | 4 m | 4 m |
| Tank length | 48 m (12 wavelengths) | 48 m (12 wavelengths) |
| Run duration | 192 s | 192 s |

There is no unique maximum height for "linear." Here the explicit criterion is
**ka <= 0.05**, selecting the largest H for the retained wavelength under that
criterion. In deep water the leading second harmonic is approximately ka/2 =
2.5% of the fundamental, while nonlinear phase-speed corrections are O(ka^2).
This increases the old linear height by about 2.49 times without making it a
steep-wave case. A stricter linearity tolerance would require a smaller H.
Residual weak nonlinearity can still produce an accuracy floor relative to an
Airy target, especially over a long propagation distance.

Depth is reduced from 18 m to one wavelength: kd = 2*pi, still effectively
deep water. The upper boundary is z=2 m, well above the waves. Generation
occupies x=0--4 m, and an 8 m numerical beach occupies x=40--48 m with length
factor 1. Both families start from a flat interface and ramp forcing over 3.2 s.

For affordability, these are **quasi-2-D extrusions**, not spanwise-convergence
studies. The periodic y extent contains two base cells (eight at the finest
level). Its physical width is 2*dx_base and changes with mesh resolution, while
dx=dy at every level. Two static factor-two refinement levels cover the wave
surface: level 1 is requested over -4H <= z <= 4H and level 2 over
-2H <= z <= 2H, across the full tank and span. Refinement boxes may expand to
cell/block boundaries. Full-volume plot and checkpoint output is disabled;
enable it explicitly when needed. The finest linear and highest-aspect-ratio
runs remain the most expensive.

## Three sweeps, for each wave model

All quoted cell sizes are on the **finest surface mesh**. The baseline uses
isotropic cells, approximately five cells per height, and nominal CFL 0.5.

### 1. Mesh resolution

| Target H/dz | Linear dx=dy=dz [m] | Actual H/dz | Stokes5 dx=dy=dz [m] | Actual H/dz |
|---:|---:|---:|---:|---:|
| 1 | 0.0625 | 1.019 | 0.375 | 0.956 |
| 2.5 | 0.025 | 2.546 | 0.15 | 2.390 |
| 5 | 0.0125 | 5.093 | 0.075 | 4.780 |
| 10 | 0.00625 | 10.186 | 0.0375 | 9.560 |
| 20 | 0.003125 | 20.372 | 0.01875 | 19.120 |

Grid-aligned spacings keep all resolutions within 5% of the requested targets.
The fixed timestep scales with cell size to preserve expected CFL 0.5.

### 2. Aspect ratio

Use the five-cells-per-height baseline's dx=dy and vary **dx/dz = 1, 2, 4, 8**.
Only dz changes; the lateral spacings and physical span remain fixed. Reduce dt
with dz to preserve expected CFL 0.5. H/dz consequently increases with aspect
ratio; it is not held at five.

### 3. CFL

Use the same isotropic baseline mesh and vary nominal CFL through
**0.1, 0.2, 0.5, 0.8, 1.0**, changing only dt and the sampling interval in steps.
The physical sampling cadence remains close to 0.04 s.

The shared baseline `mesh_cph5` belongs to all three sweeps: there are
**12 unique runs per family, 24 total**, not three copies of each baseline.
Every directory includes `input.inp`, its own `static_box.txt`, and `case.json`.
[The generated table](input_files/sweep_summary.csv) lists every actual spacing,
timestep, CFL target, and sampling cadence.

### CFL definition and limitations

For inviscid VOF with force-CFL disabled, Kynema-SGF's interface convention is

```text
C_est = 2 dt max_phase (|u_surface|/dx + |w_surface|/dz).
```

The generator selects dt = C_requested / estimated_rate. Linear velocities
are evaluated on the Airy surface with finite-depth vertical dependence.
Stokes5 planning velocities use the deep-water limit of Fenton's fifth-order,
H-based expansion evaluated on its surface. At kd=2*pi this is an appropriate
planning approximation; **accuracy references do not use that approximation**.
Phase maxima are sampled over 8192 points. The solver also checks bulk velocity
Courant numbers; its doubled interface rate is expected to dominate here.
This is not phase-speed Courant number c*dt/dx.

These are **expected, not guaranteed measured CFLs**. Air motion, interface-cell
extrapolation, reflection, startup and spatial discretization can change the
actual maximum. `time.fixed_dt` and `time.initial_dt` are equal; `time.cfl` records
the target but does not adapt a fixed timestep. Check reported CFL after startup
and calibrate dt if exact achieved-CFL matching is required. In particular,
the CFL=1 cases sit near the VOF limit and may fail if the estimate is low.
Do not silently reduce dt in a failed run and keep the original metadata.

The conventions and target-field sampling options were checked against
[Kynema-SGF's CFL calculation](https://github.com/kynema/kynema-sgf/blob/main/src/incflo_compute_dt.cpp),
[FreeSurfaceSampler](https://github.com/kynema/kynema-sgf/blob/main/src/utilities/sampling/FreeSurfaceSampler.cpp),
and [wave forcing documentation](https://github.com/kynema/kynema-sgf/blob/main/docs/sphinx/user/inputs_ocean_waves.rst).
The Stokes timestep estimate follows the deep-water limits of Fenton (1985),
with the corrected coefficients referenced in the solver's
[Stokes implementation](https://github.com/kynema/kynema-sgf/blob/main/src/ocean_waves/relaxation_zones/stokes_waves_K.H).

## Generate, run, analyze

Commands below assume the current directory is `Kynema_SGF_cases/regular_waves`.
Generation uses Python 3; analysis and tests additionally require NumPy.

```sh
python3 generate_cases.py

# Copy a self-contained case to its run directory.
mkdir -p runs/linear/mesh_cph5
cp input_files/linear/mesh_cph5/input.inp \
   input_files/linear/mesh_cph5/static_box.txt \
   input_files/linear/mesh_cph5/case.json runs/linear/mesh_cph5/

# Replace kynema_sgf with the executable / MPI launch used on your machine.
(cd runs/linear/mesh_cph5 && kynema_sgf input.inp > run.log 2>&1)
python3 analyze_case.py runs/linear/mesh_cph5

# Optional harmonic/phase figures require Matplotlib.
python3 analyze_case.py runs/linear/mesh_cph5 --plot
```

Repeat for the other generated directories. Always launch from the case's run
directory so the refinement file resolves and outputs do not collide. Keep
`case.json` with the run; if changing input values, update metadata consistently.
Regeneration overwrites generated inputs and metadata, not simulation results.

### Shared sampling

Two synchronized `FreeSurfaceSampler` groups record the same **561 points along
x=8--36 m at y=0** (0.05 m spacing, 80 samples per wavelength):

* `sampling_fs_domain` uses `phase_field = vof` for the computed surface.
* `sampling_fs_target` uses `phase_field = ow_vof` for the prescribed target.

Setting sampled `fields = ow_vof` alone would **not** locate the target surface;
the explicit `phase_field` is essential. No beach-interpolation option is used,
and all sampling points lie outside the forcing and damping regions.
The physical output interval is approximately 0.04 s, rounded to at least one
fixed timestep. Coarsest runs sample every step, giving intervals up to 0.086 s.
Sampler resolution is not CFD resolution: oversampling a coarse mesh does not
recover unresolved harmonics or interface shape.

ASCII profiles have a five-line header; the reader sorts unordered rows and
uses the x,y,z columns. The configured upper x endpoint is dropped, leaving
560 spatial samples for a periodic-length 28 m FFT. Time is inferred as
`step * dt_s` from filenames and metadata; adaptive timesteps or time-offset
restarts are not supported. Modern solver `_info.txt` sidecars are ignored.

### Accuracy and reflection from every run

The default window is **112 <= t < 176 s**, approximately 40 wave periods,
with 0.015625 Hz nominal frequency resolution and 2*pi/28 rad/m wavenumber
resolution. Startup is excluded: even a beach reflection returning to the
upstream end of the sampling line has had time to traverse the tank.
This is a defensible starting window, not a proof of stationarity. Compare
two approximately 20-period subwindows:

```sh
python3 analyze_case.py runs/linear/mesh_cph5 --start 112 --end 144 \
  --output runs/linear/mesh_cph5/analysis/early
python3 analyze_case.py runs/linear/mesh_cph5 --start 144 --end 176 \
  --output runs/linear/mesh_cph5/analysis/late
```

Extend the run if those results disagree materially.

`analyze_case.py` writes `analysis/results.json` and `analysis/results.npz`:

* **Accuracy:** space-time L1, L2, Linf and bias against the synchronized target,
  plus norms normalized by H and spatially resolved RMSE. No phase alignment or
  fitted frequency is used, so phase error remains in the norms. This is a
  **same-grid prescribed-target comparison**, including propagation/reflection
  error but not isolating the target reconstruction's own interpolation error
  against a continuum analytical surface.
* **Reflection:** remove the temporal mean at each x, apply periodic Hann
  windows in space/time, perform the 2-D transform with kernel
  `exp(-i k x + i omega t)`, and integrate directional elevation variance.
  At positive omega, k>0 is incident (+x) and k<0 reflected (-x).
  `Cr = sqrt(E_reflected/E_incident)` is the amplitude coefficient; its square
  is the variance ratio. The **0.4--0.85 Hz** band includes the fundamental and
  excludes Stokes bound harmonics. Zero-k and Nyquist modes are not directional.
* **Leakage diagnostic:** the same decomposition is applied to the prescribed
  target. Its apparent reflection helps identify sampling/window leakage, but
  is not subtracted from the simulation's coefficient.
* **Harmonics and phase:** jointly fit a mean and harmonics 1--5 to computed
  and target histories at every x. Report harmonic amplitudes, amplitude errors,
  wrapped phase errors, equivalent time lags and fit residuals, as described
  below.

The NPZ contains sorted `k_rad_per_m`, positive `omega_rad_per_s`,
`variance_m2`, `target_variance_m2`, `x_m`, and `rmse_by_x_m`. Spectral arrays have
shape `(number of omega bins, number of k bins)` and contain **variance per bin**,
not density; sum bins to integrate. Missing snapshots or target histories,
nonuniform/changing coordinates, nonfinite values, and uncovered windows raise
errors rather than returning an accuracy fallback.

### Harmonic amplitudes and phase errors

The target alone determines the fundamental frequency: locate its spectral peak
in the fundamental band, then refine by joint harmonic least squares. This
avoids locking to the nearest FFT bin, particularly important for short records.
For speed, frequency estimation uses up to 16 spatial probes; the final fits
use **every** x point. The computed surface is fitted at that same frequency,
never at an independently optimized frequency or an aligned phase.
Consequently, a simulation frequency error remains as phase error and a
nonzero fit residual.

The model is

```text
eta(x,t) = mean(x) + Re[sum_h C_h(x) exp(-i h omega (t-t0))].
A_h(x) = |C_h(x)|
phase_error_h(x) = arg[C_computed,h(x) conj(C_target,h(x))]
```

Fits are unwindowed and simultaneous, including the mean, so non-integer-period
records do not leak one fitted harmonic into another for signals matching the
model. This does not remove noise, startup, free waves or nonlinear modulation;
the residual measures what a stationary harmonic model fails to explain.
There is no directional separation within these local harmonic coefficients:
their phases contain both incident and reflected contributions.

**Positive phase error means a temporal lag at a fixed x**, with the stated
negative-time-exponential convention. Errors are wrapped to [-180,180) degrees;
reported circular means use target-amplitude-squared spatial weights, not an
arithmetic mean of angles. Equivalent lag is phase_error/(h*omega), modulo that
harmonic's period; it is not an unwrapped cumulative travel-time error. Circular
resultant length reports phase coherence across x, while RMS and maximum
absolute phase errors show variation hidden by a mean.

Phase is undefined where either fitted amplitude is below **0.001 times the
spatial RMS target fundamental amplitude**. Such points are masked; unresolved
phase/error ratios are JSON `null` with an explicit status. NPZ phase arrays use
NaN and an accompanying validity mask. For example, an ideal linear target has
no meaningful phase reference for harmonics 2--5, even if the computed signal
contains those harmonics. Their measured amplitudes are still reported.
The floor is an amplitude screen, not a statistical confidence interval.
Compare small harmonics with fit residuals and subwindow sensitivity before
interpreting their phases.

JSON `harmonic_analysis` contains the fundamental frequency, phase convention,
threshold, residual norms, and a `harmonics` list. Each entry reports spatial RMS
amplitude (not RMS elevation), amplitude normalized by the target fundamental,
relative amplitude error where resolved, circular-mean/RMS/maximum phase error,
valid point count, phase resultant length and signed equivalent lag.

The NPZ additionally contains `harmonic_orders`, `harmonic_frequencies_hz`,
`harmonic_signal_coefficient_m`, `harmonic_target_coefficient_m`,
`harmonic_signal_amplitude_m`, `harmonic_target_amplitude_m`,
`harmonic_phase_error_rad` and `harmonic_phase_valid` with shape
`(number of harmonics, number of x points)`, plus the per-x fit residuals.
Orders/frequencies are 1-D. Complex coefficients retain the reference phase
at `time_origin_s`, shared between computed and target fits.

Use `--harmonics N` or `--phase-floor-ratio R` to change those choices. The
highest searched harmonic must stay below sampling Nyquist, and insufficient
samples, an ill-conditioned fit, a frequency-band-edge solution, or an absent
target fundamental raise explicit errors. `--plot` adds
`analysis/results_harmonics.png`: amplitude/phase by order and local phase error
against x. Matplotlib is optional unless plots are requested.

## Historical data

The second-order inputs, archived records and derived outputs have been removed.
Remaining `data/reflection` and `data/convergence` records are historical;
they were **not** produced by the new inputs (notably H, depth, domain and dt
differ). Existing Julia scripts and standalone k-omega CLI defaults refer only
to those archived runs. Use `analyze_case.py` for all new cases. Historical
records without target histories cannot provide the new paired analysis.

### Reproducible sample visualizations from existing data

The archived fifth-order **convergence** runs include a target gauge history at
x=10 m, allowing genuine harmonic and phase comparisons for all five legacy
resolutions. This reference is shared across the archive and is **not assumed
same-grid** for every legacy run. Unlike new cases, these runs retain historical
geometry and a common fixed dt rather than matched CFL.

```sh
python3 plot_archived_harmonics.py
```

This requires NumPy and Matplotlib and writes the following to
`plotting_outputs/`:

* `archived_harmonics_phase.png`: harmonic amplitudes, phase errors and raw
  computed-minus-target histories;
* `archived_harmonics_convergence.png`: time traces, amplitude and phase errors
  against mesh resolution, and the finest run's harmonic-fit residual;
* `archived_harmonics.json`, `.csv` and `.npz`: metrics, harmonic coefficients,
  masks and sampled histories behind the figures.

The default window **10 <= t < 18 s** contains only about 5.2 periods; the fitted
target frequency is approximately **0.650018 Hz**. Fundamental diagnostics are:

| Archived H/dx | Fundamental amplitude error | Phase error (positive = lag) |
|---:|---:|---:|
| 1.43 | -37.15% | -33.62 deg |
| 2.87 | -15.98% | -19.74 deg |
| 5.74 | -7.08% | -10.58 deg |
| 11.47 | -2.91% | -5.33 deg |
| 22.94 | -1.18% | -2.66 deg |

These are preliminary **archived-data** diagnostics, not results from the new
24-case sweeps and not uncertainty bounds. Negative phase indicates a lead.
Higher harmonics, especially 4--5, show substantially less uniform improvement
and are more sensitive to residuals and the short fitting window. No analytic
target or phase error is invented for archived linear reflection profiles that
lack a paired stored target; the new linear cases explicitly collect that data.
Use `--start`, `--end`, `--data` and `--output` to reproduce other archived windows
or save a comparison separately.

## Tests

```sh
PYTHONPATH=reflection python3 -m unittest \
  test_regular_wave_cases test_reflection_coeff_komega test_harmonic_analysis
```

Tests check actual rendered mesh spacings and cell-height targets, lateral
isotropy, aspect ratios, CFL/dt scaling including the interface factor of two,
complete sweep membership with shared baselines, reproducible generation,
paired surface/target analysis, known synthetic reflection ratios, off-bin joint
harmonic recovery, signed phase/lag, circular wrapping, weak-component masking,
retained signal-frequency error, harmonic exclusion, leakage, malformed
histories, and persistent output files. The plot smoke test runs when Matplotlib
is installed.
Solver execution is a separate validation step; these tests do not establish
measured CFL, stability or stationarity.
