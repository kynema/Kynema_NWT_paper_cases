# Irregular periodic ocean post-processing

## Instantaneous kx-ky comparison

From the repository root, using the NumPy and Matplotlib dependencies declared in
[`tools/ocean/pyproject.toml`](../../tools/ocean/pyproject.toml):

```sh
python3 Kynema_SGF_cases/irregular_ocean_periodic/instantaneous_spectra.py
```

The script also works from other working directories. It writes
[`plotting_outputs/instantaneous_2D_spectra.png`](plotting_outputs/instantaneous_2D_spectra.png).
The same comparison is included in the existing
[`spectral_analysis.py`](spectral_analysis.py) marimo notebook; existing 1D spectra
and energy plots are unchanged.

Rows show HOS-Ocean, all five SGF configurations, and the target solution.
Columns show the available snapshots in chronological order (currently initial
and approximately 1000 s). Titles report the actual file-header times:
HOS ends at 999.6 s, while SGF/target end at 999.989 or 1000 s. The script rejects
missing snapshot sets or chronological matches differing by more than 1 s.

All panels share physical wavenumber axes, equal aspect ratio, and a single
logarithmic color scale. The displayed quantity is surface-mode amplitude in
metres, **not power or spectral density**. HOS amplitudes and coordinates are
read directly from the split `a_3d_time*.txt` files, with shared coordinates
from `a_3d_time0.txt`. SGF/target elevations are transformed using a 2D FFT on
their uniformly spaced periodic sampling grids. Positive-kx coefficients are
doubled except at the even-grid Nyquist column, matching the
[HOS Fourier normalization](https://github.com/LHEEA/HOS-ocean/blob/master/sources/utilities/fourier_r2c_FFTW3.f90).
No window, mean subtraction, or per-panel rescaling is applied.
Values at or below 1e-6 m (including zero) use the lowest color; the upper
limit is the maximum across every panel. Signed ky is retained; NumPy places
the even-grid ky Nyquist mode at the negative edge, whereas HOS uses the positive edge.

Run the numerical and plotting checks with:

```sh
MPLBACKEND=Agg python3 -m unittest discover \
  -s Kynema_SGF_cases/irregular_ocean_periodic -p 'test_instantaneous_spectra.py'
```
