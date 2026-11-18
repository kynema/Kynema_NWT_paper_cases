# Four-jacket Morison load check

This is an executable **rigid, normal-strip Morison benchmark**, not an OpenFAST
structural simulation or a replacement for full HydroDyn. It uses the published
jacket diameters, as requested, and does not modify the CFD mesh.

For actual standalone HydroDyn input decks, including the selected flooded-leg
baseline, hydrostatics and axial/end loads, use the
[version-pinned HydroDyn setup](../hydrodyn/README.md). Its collected loads can
also be analyzed by this directory's comparison command.

## Geometry and provenance

The source is Mroczek, Arwade and Lackner (2023),
[Wind Energy Science 8, 807-817](https://doi.org/10.5194/wes-8-807-2023).
The paper links to
[`mmrocze2/IEA-15-240-RWT`, commit
`da77a393569b6c5ddf7aa424e5b0a92d37eba833`](https://github.com/mmrocze2/IEA-15-240-RWT/tree/da77a393569b6c5ddf7aa424e5b0a92d37eba833/OpenFAST/IEA-15-240-RWT-Jacket_60m_water_depth).
An unchanged copy of its HydroDyn file and upstream Apache-2.0 license are in
[reference/](./reference/). The parser checks the source SHA-256 before extracting
27 joints, 57 members and 9 diameter sets. Source SHA-256:
`a017b4d7711e3673e88f520dfe4f61a88e2086af10169607736319c68885549d`.

**Rebuild the mesh before using this as an accuracy validation.**
The archived CFD wall is consistent with approximately **twice the published
tube diameters**. For example, the published pile-stub diameter is 3 m, whereas
the current mesh has approximately 3 m **radius**. The saved
[geometry audit](./geometry_audit.json) records the mesh checksum and residuals:
median distance from wall nodes to published cylinder surfaces is 0.517 m,
versus 0.012 m using doubled diameters. Joint intersections and caps prevent
this residual from being zero; the audit is not proof of complete CAD identity.

For remeshing, use **radius = HydroDyn `PropD` / 2**. Source coordinates have
SWL at z=0. Mesh-local coordinates are source coordinates plus `(0,0,60)`;
clip buried pile ends at local z=0. Retain the four existing translations
`(-390,-390,-60)`, `(-390,390,-60)`, `(390,-390,-60)`, `(390,390,-60)`,
with no rotations. Thus world coordinates are source coordinates plus each
jacket's **horizontal** displacement; do not subtract another 60 m.

The vendored HydroDyn file is a source model, **not a generated modern
standalone-HydroDyn deck**. The paper used modified OpenFAST 3.2.1
([openfast/openfast#1008](https://github.com/OpenFAST/openfast/pull/1008)).
Its flooding-group line also appears inconsistent (count 1 with 15 member
IDs). No automatic repair or version conversion is attempted here.

## Model and interpretation

For each fixed circular member, with axis unit vector `s` and projector
`P = I - s s^T`, the force per unit length is

```text
uN = P u
aN = P (partial u / partial t)
f  = 0.5 rho Cd D |uN| uN + rho Cm (pi D^2 / 4) aN
M  = integral ((x - base_origin) cross f) dl
```

Defaults are `rho=1025 kg/m^3`, `Cd=1.05`, and **`Cm=Cp+Ca=1.0+0.2=1.2`**.
These are an explicit coefficient baseline from the source, not a calibration
to smooth CFD tubes. Repeat with justified alternative `--cd`/`--cm` values;
do not independently tune each jacket.

Members are subdivided to at most 0.5 m by default, respecting any smaller
source `MDivSize`. End elevations define a linearly interpolated surface.
Strips are clipped to that surface and the z=-60 m seabed; loads are evaluated
at the resulting wet midpoints. Refine `--spacing` to check strip-integration
convergence. This is centerline wetting, not partially submerged cross-section
integration or a breaking-wave/slamming model.

Velocity and acceleration must be **undisturbed incident** fields with the same
directional phases, physical clock and current as CFD. `hos` uses the full
HOS velocity field, but its inertia term uses a fixed-position **Eulerian**
time derivative, a conventional leading-order engineering approximation.
It does **not** implement nonlinear pressure-gradient/Froude-Krylov integration
or convective acceleration. It never differentiates velocity along the moving
wet midpoint; that would introduce a spurious acceleration.

Excluded: hydrostatics, axial/end loads, flooding/internal-water effects, marine
growth, joint-overlap corrections, diffraction/shielding, structural motion,
wind, soil and global damping. In particular:

- The paper's growth/flooding and structural damping are not silently carried
  into the bare, fixed CFD comparison.
- This is **not total six-component HydroDyn equivalence**. All six normal-strip
  load components are saved, but total CFD Fz cannot be validated with this model.
  Time-varying buoyancy/end loads can also affect overturning moments.
- Demeaning removes a constant offset only; it does not fix those missing terms.
  Inspect both raw and demeaned errors rather than calling a mean removal a
  buoyancy correction.
- For a full force-budget check, separately validate a version-compatible
  HydroDyn model with the intended flooding/caps, hydrostatics and axial terms.
  Do not simply run the archived source file with an arbitrary OpenFAST version.

## Setup

Commands below run from the `offshore_array` directory. Python >=3.9 and NumPy
>=2.0 are sufficient for integration/tests. Matplotlib is needed only for
`--plot`; `netCDF4` only for `check-mesh`. Install these in your normal scientific
Python environment as needed. No OpenFAST executable is required by these scripts.

```bash
python3 morison/run.py check-mesh input_files/jacket_mesh_origin.exo \
  --output morison/geometry_audit.json
# Exit 2 means a geometry mismatch, as expected for the current archived mesh.
# Re-run after remeshing; a radius-check pass is not a full CAD validation.

python3 morison/run.py prepare --spacing 0.5 --output morison_outputs/prepared
```

This exports a geometry NPZ, provenance JSON and `strip_endpoints.csv`, with
world coordinates and jacket/member IDs. Geometry hashes prevent inadvertently
mixing kinematics from different strip spacings or coordinate transforms.

### 1. Archive and register the actual incident realization

Required real inputs are **not included** in this case directory:

1. The demonstration's own `modes_HOS_SWENSE.dat` and HOS input/log files,
   currently referenced by the SGF input under
   `/scratch/mkuhn/HOS_Ocean_precursor_halflength_current/Results/`.
   Do not substitute the separate periodic-ocean example's modes or a
   random-phase JONSWAP realization.
2. The rebuilt-mesh CFD `jacket1_forces.dat` through `jacket4_forces.dat`,
   with original headers, physical timestamps, moment origins, force/pressure
   sign and pressure convention. Preserve logs, restart provenance and mesh.
3. Paired incident-field registration records from the updated SGF inputs.

Both [fresh](../input_files/offshore_array_sgf.inp) and
[restart](../input_files/offshore_array_sgf_rst.inp) inputs now save 12 point
probes every five SGF steps: three depths (-2, -15, -40 m) at each jacket center.
They record `ow_velocity`, `ow_levelset`, CFD `velocity`, and `vof` in ASCII.
The [point file](../input_files/morison_registration_points.txt) is relative to
the case working directory, like the existing mesh paths. Keep actual sample
times: the coupled solver uses adaptive stepping.

Use the **`ow_*` target fields**, not structure-disturbed CFD velocity, to verify
HOS/CFD spatial and temporal registration. Compare multiple locations and wave
cycles; use deeper probes when the -2 m probe is dry. Check the level-set sign
when recovering elevations. Retain the CFD fields to diagnose incident-field
versus local-flow differences, not to drive the incident Morison calculation.

The configured HOS offsets are -750 m in x/y, suggesting
`HOS_xy = CFD_world_xy + (750,750)`; **verify this with probes for your executable**.
Also check HOS initialization time, wave direction, velocity units, and any
current included in the exported modes. A directory name containing `current`
does not establish the velocity convention. No current is silently added here.
Check logs for HOS end-of-file/recycling warnings: a recycled precursor is
not the same physical time record. The adapter rejects times beyond its record.

### 2. Compute Morison histories from HOS

The optional adapter targets the documented GNU-Fortran module interface of
[LHEEA/Grid2Grid at commit
`f5f997a35b940dd6c479b88019c2087d98d6144d`](https://github.com/LHEEA/Grid2Grid/tree/f5f997a35b940dd6c479b88019c2087d98d6144d).
Build its shared library using the upstream instructions (Fortran/FFTW and its
other build dependencies). The target ABI is documented in its
[interface example](https://github.com/LHEEA/Grid2Grid/blob/f5f997a35b940dd6c479b88019c2087d98d6144d/doc/interface.tex).
Other versions/compilers are not assumed ABI-compatible. Grid2Grid is an
external dependency, not vendored here.

Before asserting registration, the same provider can be used to inspect
incident velocities directly:

```python
import sys
import numpy as np
sys.path.insert(0, "morison")
from hos import Grid2Grid

hos = Grid2Grid("/absolute/path/libGrid2Grid.so", "/absolute/path/modes_HOS_SWENSE.dat")
points = np.loadtxt("input_files/morison_registration_points.txt", skiprows=1)
# Example candidate mapping only; compare at the recorded SGF sample times.
time = 240.0
print(hos.elevation(points[:, :2] + [750, 750], time))
print(hos.velocity(points + [750, 750, 0], time))
```

After confirming registration, run, for example:

```bash
python3 morison/run.py hos morison_outputs/prepared/geometry.npz \
  --library /absolute/path/libGrid2Grid.so \
  --modes /absolute/path/modes_HOS_SWENSE.dat \
  --hos-xy-shift 750 750 --hos-time-shift 0 --registration-confirmed \
  --start 240 --stop 3600 --dt 0.15 --derivative-step 0.02 \
  --output morison_outputs/morison_loads.npz
```

The example requires modes covering at least `[239.98,3600.02]` s. Adjust the
window to the **recorded common interval**. The adapter requires 60 m HOS depth.
It queries surface elevations at strip ends, then velocities at wet midpoints,
and computes centered differences at **fixed coordinates**. Dry strips need
no velocity queries. Loads are integrated frame-by-frame without storing a
huge time-by-point kinematics archive.

Repeat with half `--derivative-step`, finer `--nz`, smaller `--dt` and smaller
strip spacing until relevant RMS/peak loads and phases are stable. Grid2Grid
and Waves2AMR use different interpolators; quantify their incident-field
differences with the probes before attributing discrepancies to CFD.

### Alternative: use another incident-field exporter

`integrate` accepts a non-pickled NumPy NPZ with:

| Key | Shape / meaning |
|---|---|
| `time` | `(nt,)`, strictly increasing physical seconds |
| `eta` | `(nt,nstrip,2)`, endpoint surface elevations, world z |
| `velocity` | `(nt,nstrip,3)`, incident velocity at each **wet midpoint**, m/s |
| `acceleration` | same shape, fixed-position partial du/dt, m/s^2 |
| `acceleration_convention` | scalar string `"eulerian"` |
| `geometry_sha256` | scalar string from prepared `geometry.json` |
| `source` | scalar string identifying actual source file/run and mapping |

Use `model.wet_geometry(geometry, eta_at_time)` to get the required sampling
positions after determining endpoint elevations. Dry-strip velocities and
accelerations may be finite zeros. Do not sample at the original dry midpoint
and pass those values as wet-midpoint kinematics. All arrays must be finite.
This interface also supports incident kinematics exported from a matching
**no-structure** SGF run; native disturbed CFD velocities are not appropriate.

```bash
python3 morison/run.py integrate morison_outputs/prepared/geometry.npz \
  incident_kinematics.npz --output morison_outputs/morison_loads.npz
```

### 3. Compare rebuilt-mesh CFD loads

```bash
python3 morison/run.py compare morison_outputs/morison_loads.npz \
  --cfd 1 jacket1_forces.dat --cfd 2 jacket2_forces.dat \
  --cfd 3 jacket3_forces.dat --cfd 4 jacket4_forces.dat \
  --cfd-format nalu10 --force-on-structure-sign 1 \
  --start 240 --dt 0.15 --plot --output morison_outputs/comparison
```

`nalu10` is an **explicit format selection**, following the legacy Julia
script's convention: time; three pressure-force columns; three viscous-force
columns; three total-moment columns. Pressure and viscous forces are summed.
The reader tolerates leading text headers but requires exactly ten numeric
columns. Check real headers before selecting this format and sign. It does
not guess different output layouts or automatically reconstruct pressure.
Alternatively supply CSV with exactly:
`time,Fx,Fy,Fz,Mx,My,Mz`, in N and N m, force **on the structure**.

All CFD moments must be about the configured jacket base centers, in world
axes. The array resultant shifts moments to world `(0,0,0)` **before** summing.
Duplicate/backward times are errors: explicitly stitch and document restart
branches first; no timestamps are jittered, silently sorted or deduplicated.
There is no fitted time lag, phase correction or amplitude rescaling.

Outputs include raw/mean-removed RMSE, means, standard deviations, extrema,
correlations, up to three separated largest-horizontal-load event windows per
body, and four-jacket plus array histories. The NPZ also contains one-sided
Hann/Welch PSDs, coherence and wrapped cross-phase for all components.
Segments contain up to 4096 samples (614.4 s at dt=0.15 s), with 50% overlap;
shorter records use at most half the available samples. The report records
segment length/count so low-frequency resolution and averaging are explicit.
Positive phase means CFD leads Morison; ignore phase in low-energy or
low-coherence bins. Spectral columns follow body-major order, six channels
`Fx,Fy,Fz,Mx,My,Mz` per body; the fifth body is the array at world origin.
Undefined zero-power phase/coherence are NaN in NPZ, not fabricated zeros.
Optional figures show time histories, force spectra and jacket-1 event windows.

Both records are linearly interpolated onto the common interval at `--dt`;
there is no extrapolation or anti-alias filtering. Preserve native histories
and choose an adequately resolved common cadence instead of strongly
downsampling broadband/slamming loads. Statistics describe this common grid.
Start with the existing 240 s startup exclusion, then verify the actual
wave/current field has settled. Use the same physical interval for every jacket.

## Verification and current limits

```bash
python3 -W error::RuntimeWarning -m unittest discover -s morison -p 'test_*.py' -v
python3 morison/run.py demo --output morison_outputs/synthetic_demo --plot
```

The demo is **synthetic only**: analytic directional linear-wave kinematics
drive the model, and the fabricated "CFD" files are exactly 1.1 times its loads.
It checks the complete file/comparison/plot workflow, not CFD agreement.
Tests also cover analytic current/inertia loads, inclined members, partial
wetting, moment transport, fixed-position differentiation, restart rejection
and the cross-phase sign.

Validation performed: all 10 focused tests passed with numerical RuntimeWarnings
treated as errors; the synthetic reports and figures were generated and checked.
The current mesh audit returned the expected mismatch (exit 2).
No actual demonstration HOS modes, CFD force histories, compatible Grid2Grid
library were available for a real-load run. The separate standalone HydroDyn
setup has since been tested with a locally built v5.0.0 executable and synthetic
fields; it has not been run with the actual demonstration's incident data.
The Grid2Grid ABI adapter and changed CFD sampling inputs still require an
on-system executable smoke test; analytic provider tests are not a substitute.
