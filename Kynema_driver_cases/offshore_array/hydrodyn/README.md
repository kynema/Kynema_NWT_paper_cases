# Standalone HydroDyn jacket inputs

These are actual input decks for **OpenFAST v5.0.0's `hydrodyn_driver`**,
not the older paper-linked file passed to a newer executable unchanged.
The [input_files/](./input_files/) directory contains:

- [HydroDyn.dat](./input_files/HydroDyn.dat): shared published-diameter jacket.
- Four drivers, [jacket1.dvr](./input_files/jacket1.dvr) through
  [jacket4.dvr](./input_files/jacket4.dvr), each referencing its own SeaState file
  and incident-data root.
- [smoke_regular.dvr](./input_files/smoke_regular.dvr): a runnable synthetic
  regular-wave check needing no external wave data.
- [case.json](./input_files/case.json): exact grid, clock, origins and provenance.
- The upstream jacket model's Apache-2.0 [LICENSE](./input_files/LICENSE).

**Selected baseline: flooded legs/pile stubs, members 1-15**, with
`FillFSLoc=0` and `FillDens=1025 kg/m^3`. Tubes are bare, without marine growth.
The source's inconsistent filled-group count is replaced by an explicit new
15-member group in the converted deck; the archived source is not edited.

## Physical scope

Geometry comes from the paper-linked
[`mmrocze2/IEA-15-240-RWT` commit
`da77a393569b6c5ddf7aa424e5b0a92d37eba833`](https://github.com/mmrocze2/IEA-15-240-RWT/tree/da77a393569b6c5ddf7aa424e5b0a92d37eba833/OpenFAST/IEA-15-240-RWT-Jacket_60m_water_depth).
See the [existing geometry audit](../morison/geometry_audit.json):
the CFD mesh still needs rebuilding with the **published diameters**.
The 0.1 m buried pile ends are clipped to z=-60 m; other joint coordinates,
diameters, wall thicknesses and connectivity are preserved.

The jacket is fixed: prescribed motions, velocities and accelerations are all
zero. Water depth is 60 m, density 1025 kg/m^3, and gravity 9.81 m/s^2.
Source normal coefficients are retained: `Cd=1.05`, `Ca=0.2`, `Cp=1`
(`Cm=1.2`), along with its distributed axial coefficients and joint `AxCp=1`.
Hydrostatic member loads are enabled (`MHstLMod=1`, `HstMod=1`), with added mass
to the instantaneous surface (`AMMod=1`). No WAMIT files, structural model,
turbine, wind, soil model, global damping or OpenFAST `.fst` file are required.

This includes HydroDyn's hydrostatic and axial/end contributions, unlike the
[normal-strip Python benchmark](../morison/README.md). It is still an engineering
member model: source joint overlaps are not removed, diffraction/shielding and
inter-jacket interaction are absent, and it does not resolve fluid inside the
jacket. HydroDyn's SWL internal fill model is not a resolved free-flooding
calculation. Confirm the rebuilt CFD model's flooding, caps/openings and joint
intersections before interpreting total Fz or buoyancy differences as CFD error.
The paper's marine growth and structural damping are deliberately not imported.

`WaveMod=6` preserves the externally supplied directional realization rather
than drawing a new random sea or reconstructing it from one elevation gauge.
SeaState's grid extends from the seabed to SWL. `WaveStMod=1` uses vertical
stretching: below-SWL kinematics are unchanged, and above-SWL wet nodes use
SWL values. This is **not exact nonlinear HOS kinematics above SWL**; quantify
that approximation for the intended sea state. Current must already be present
in the supplied fields; SeaState current and second-order synthesis are off.

## Run now: executable smoke test

Commands assume the `offshore_array` working directory. Use a
**v5.0.0 standalone driver**, not the OpenFAST turbine executable:

```bash
/path/to/hydrodyn_driver hydrodyn/input_files/smoke_regular.dvr
```

This runs 36.75 s of synthetic H=2 m, Tp=12.25 s regular waves on one fixed
jacket. It writes `smoke_regular.out`, echo files and module summaries beside
the inputs. This is an installation/input-format test, **not demonstration data**.

The decks were checked against
[OpenFAST v5.0.0](https://github.com/OpenFAST/openfast/tree/v5.0.0) and its pinned
[regression inputs](https://github.com/OpenFAST/r-test/tree/dd5feaaaa500ba7283140107806300d551cff0a7/modules/hydrodyn/hd_MHstLMod2).
The standalone target can be built with:

```bash
cmake -S /path/to/openfast-v5.0.0 -B /path/to/build \
  -DBUILD_TESTING=OFF -DOPENMP=OFF -DCMAKE_BUILD_TYPE=Release
cmake --build /path/to/build --target hydrodyn_driver --parallel 2
```

A Fortran compiler, CMake and BLAS/LAPACK are required. On the local macOS
validation build, `-DBLA_VENDOR=Apple` selected the existing Accelerate libraries.
The executable is under `build/modules/hydrodyn/`. Source-archive builds may
print `HEAD-HASH-NOTFOUND` because there is no Git checkout; use the pinned
release source, not that banner, to identify this validation build.

## Run later: the actual four-jacket realization

Generate a working copy rather than editing templates or generating outputs
in the repository:

```bash
python3 hydrodyn/prepare.py --output hydrodyn_outputs/case
```

Defaults select **physical times 240-3675 s**, with HydroDyn local time
0-3435 s and `dt=0.15 s`. Both endpoints are included: 22,901 records.
Change `--physical-time-origin`, `--duration` and `--dt` to match the available
common HOS/CFD record; duration must be an integer multiple of dt.
All drivers have identical clocks. The generator refuses to regenerate a
directory that already contains incident kinematics.
Input-file hashes in the manifest detect inconsistent manual edits; regenerate
after changing templates or numerical settings rather than mixing modified decks
with old geometry/grid metadata.

Each jacket uses local `(x,y)` coordinates centered on that jacket and SWL z=0:

```text
world x,y = local x,y + jacket center x,y
world z   = local z
HOS x,y   = world x,y + verified HOS-domain shift
physical time = HydroDyn time + physical_time_origin
```

Do not translate geometry to +/-390 m inside HydroDyn or add another -60 m
vertical shift. Distinct local wave grids, sampled from the **same HOS
realization at the four world locations**, encode the phase differences.
Use the [incident registration probes](../input_files/morison_registration_points.txt)
and [registration guidance](../morison/README.md) to verify positions, clocks,
directions and current. Do not copy one jacket's incident files to the others.

### Information still needed

1. The actual demonstration HOS modes and input/log files, covering the selected
   physical window without precursor recycling.
2. Registered full fields on the grids in `case.json`: three velocities,
   three **Eulerian fixed-position** accelerations, dynamic pressure, and surface
   elevation. A single elevation trace is insufficient.
3. Pressure in **Pa at rho=1025**, excluding hydrostatic and atmospheric-gauge
   pressure. Check any reconstruction library's pressure units and density
   convention; HydroDyn adds its own hydrostatic loads.
4. Rebuilt-mesh CFD histories and confirmation of flooding/end conditions.

The existing Grid2Grid wrapper only exports velocity/elevation for the Python
benchmark. It does **not** yet supply this complete pressure/grid handoff.
Do not substitute zero pressure or acceleration for unavailable real data.
Obtain these fields from a validated HOS reconstruction/exporter or matching
no-structure calculation, preserving directional phases and current.

### Grid and data handoff

Defaults use a +/-25 m square local grid: 11 x 11 horizontal nodes and 31
nonuniform vertical nodes, ordered from seabed to SWL. Grid nodes are not
HydroDyn member-integration nodes. Export their world locations with:

```bash
python3 hydrodyn/wavekin.py grid hydrodyn_outputs/case \
  --jacket 1 --output hydrodyn_outputs/jacket1_grid.csv
```

Repeat for jackets 2-4. Initial grid spacing is 5 m horizontally; refine
`--nx-half`, `--ny-half`, `--nz` and `--member-ds` in new cases to check
interpolation and load convergence. It is not a demonstrated converged grid.
Default wave-array storage is **at least 2.25 GiB per jacket**, excluding other
allocations; ASCII incident files are considerably larger. Run the four jackets
sequentially unless memory permits concurrent runs.

For each jacket, prepare a directory of memory-mappable `.npy` arrays:

| File | Shape / meaning |
|---|---|
| `x.npy`, `y.npy`, `z.npy` | Exact local axes from `case.json`; z bottom-to-top |
| `time.npy` | `(nt,)`, **physical** times |
| `velocity.npy` | `(nt,nz,ny,nx,3)`, world-axis components, m/s |
| `acceleration.npy` | Same shape, fixed-position partial du/dt, m/s^2 |
| `dynamic_pressure.npy` | `(nt,nz,ny,nx)`, Pa, hydrostatic excluded |
| `elevation.npy` | `(nt,ny,nx)`, surface z relative to SWL, m |

Include `metadata.json`, identifying the actual source and verified mapping:

```json
{
  "jacket": 1,
  "registration_confirmed": true,
  "source": "Actual precursor/run, coordinate and clock mapping, reconstruction settings",
  "dynamic_pressure_units": "Pa",
  "acceleration_convention": "eulerian"
}
```

Then export and validate:

```bash
python3 hydrodyn/wavekin.py export hydrodyn_outputs/case \
  --jacket 1 --bundle /path/to/jacket1_bundle
# Repeat for jackets 2-4 with their respective bundles.
python3 hydrodyn/wavekin.py validate hydrodyn_outputs/case
```

The exporter writes each jacket's `kinematics/jacketN/waves.{Vxi,Vyi,Vzi,Axi,Ayi,Azi,DynP,Elev}`.
It memory-maps the arrays and writes one time frame at a time; it does not need
a huge compressed NPZ in memory. It checks axes, clock, finite data and
registration metadata, refuses overwrites, and marks failed exports INCOMPLETE.
Validation reads all 32 files and checks every numeric row.

Alternatively, supply native files directly, following the **v5.0.0 parser**:
each file has **exactly 13 header lines**, then x-fast rows. For each inclusive
time sample, the seven volumetric files contain `nz*ny` lines with `nx` values
per line (bottom-to-top z, then increasing y); `.Elev` contains `ny` such lines.
There are **no time/index columns**. This follows `UserWaves_Init` and
`SeaStOut_WriteWvKinFiles`; do not follow contradictory statements in some
documentation about missing headers or omitted final samples.
Nonnumeric fields can silently become zero in SeaState itself, which is why the
strict validator should be run before the executable.

### Execute and compare

After supplying and validating the real incident fields:

```bash
/path/to/hydrodyn_driver hydrodyn_outputs/case/jacket1.dvr
/path/to/hydrodyn_driver hydrodyn_outputs/case/jacket2.dvr
/path/to/hydrodyn_driver hydrodyn_outputs/case/jacket3.dvr
/path/to/hydrodyn_driver hydrodyn_outputs/case/jacket4.dvr

python3 hydrodyn/collect.py hydrodyn_outputs/case \
  --output hydrodyn_outputs/hydrodyn_loads.npz
```

`HydroFxi/Fyi/Fzi` are forces on the jacket in N. `HydroMxi/Myi/Mzi`
are in N m about the local **SWL origin**, not the seabed. The collector adds
`(0,0,60) cross F` to move moments to each CFD base, restores the physical time
origin, and checks output completion, times and channel units.
The NPZ can directly replace the model-load argument of
[`morison/run.py compare`](../morison/README.md); plots identify it as HydroDyn,
not the normal-strip Python model. Keep the same CFD format/sign and analysis
window checks described there.

## Validation performed

Built the standalone driver from the v5.0.0 source and ran:

- The regular-wave deck through all 246 output times, ending at 36.75 s.
- Four WaveMod=6 decks with separately exported synthetic constant-current
  fields; all terminated normally, with valid requested channel units. HydroDyn's
  summaries confirmed the published radii and exactly members 1-15 flooded.
- The force scaling check `Fx(jacket N)/Fx(jacket 1)=N^2` for speeds N m/s,
  plus the output collection, SWL-to-base moment shift and CFD-comparison handoff.

All 15 focused HydroDyn/Morison tests passed with numerical RuntimeWarnings
treated as errors. To repeat executable tests, set `HYDRODYN_DRIVER` to the
absolute standalone executable path:

```bash
HYDRODYN_DRIVER=/path/to/hydrodyn_driver \
  python3 -W error::RuntimeWarning -m unittest discover -s hydrodyn -p 'test_*.py' -v
```

These tests establish executable/file-format operation, **not agreement with
the real demonstration**, whose complete incident fields remain unavailable.
