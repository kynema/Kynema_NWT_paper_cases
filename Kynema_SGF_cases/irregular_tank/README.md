# Irregular tank: temporal extremes and crest statistics

`temporal_analysis.py` compares only the **HOS-NWT precursor** and the two
**Kynema-SGF** gauge records. It never loads experimental data. The legacy
`IrregularNWT.jl` time-series/PSD entrypoint also no longer requires the
experimental record or its external path.

## Run

Python 3 and NumPy are required; figures additionally require Matplotlib.
From this directory:

```sh
python3 temporal_analysis.py

# Rank by complete-wave crest-to-trough height rather than crest height.
python3 temporal_analysis.py --rank-by height --output plotting_outputs/by_wave_height

# Change the shared record and number of event windows.
python3 temporal_analysis.py --start 60 --end 1700 --events 4 \
  --output plotting_outputs/alternate_window
```

Default analysis uses **30--1776.99 s**, the shared portion of all three
records after excluding the first 30 s (about 18 peak periods). Times,
elevations and the header's **Tp = 1.6971 s** are already dimensional
**model-scale** quantities; no scale-50 experimental conversion is applied.
`--end` defaults to the earliest available record end. `--data` selects a
different directory containing the same three filenames, and `--no-plots`
produces only the numerical reports without requiring Matplotlib.

## Gauge and clock handling

* HOS uses `HOS_NWT_probes.dat`, probe **p2** (the third numeric column), matching
  the existing Julia analysis. SGF uses time and z from the two
  `amr_wind_hos_nwt_sampling_time_series*.txt` files at **x=y=0**.
* The expected coordinate mapping is HOS x=16 m to SGF x=0, using the input's
  `HOS_domain_offset_x = -16`. However, the archived HOS header refers to
  `probe_2.inp`, which is not included here. The available precursor
  `probe_1.inp` contains x=16 m but does not establish both archived positions.
  The p2 choice therefore follows prior processing; verify the missing position
  list before treating collocation as independently confirmed.
* HOS timestamps are printed with limited scientific-notation precision and
  eventually repeat. Recover the uniform clock as `t0 + index/f_out`, where
  **f_out = 141.421 Hz** comes from the file header. Every reconstructed value
  must lie within the printed timestamp's rounding interval, or analysis fails.
  Using the late-record median printed dt or adding tiny increments to ties
  would yield an incorrect clock.
* SGF's adaptive-step records also have repeated printed timestamps: 209 tied
  samples in the coarse record and 8064 in the fine record within the default
  window. **Keep all elevations and their original order for native crest
  extraction.** Only interpolation and time-weighted moments average elevations
  sharing an identical printed timestamp.
* The three histories remain on their **original physical clock**. There is no
  cross-correlation alignment, fitted shift, independent event recentering or
  dynamic time warping. Peak-time differences cannot be interpreted below the
  underlying sampling and printed-time precision.

## Largest-wave event comparisons

1. Calculate each record's time-weighted mean over the same analysis window.
2. Detect complete waves between consecutive **upcrossings of that mean**,
   using linear crossing times and native sampled crest/trough values.
   Partial waves at record edges are excluded.
3. Rank **HOS** waves by crest height above the HOS mean, or by crest-to-trough
   wave height with `--rank-by height`.
4. Choose six events by default, requiring full **+/-3 Tp** context and
   nonoverlapping windows. This separates groups visually; it is not proof of
   statistical independence.
5. Plot all native records around the **same HOS crest times**. Shading marks
   the complete HOS wave interval, and dots mark matched native maxima.

For each SGF record, a matching crest is the maximum inside the **same HOS
upcrossing interval**, not a maximum chosen from the entire context window.
Report its crest-height error, relative error, peak-time lag (positive means
later than HOS), and context-window RMSE, bias, MAE and correlation. A maximum
at the matching interval's edge is flagged because it may not be a properly
matched crest. Full-record crest statistics also show each SGF record's own
largest crests, including those not selected by the HOS event ranking.

Time-series agreement is evaluated on a shared **0.01 s linear-interpolation
grid**; `--sample-dt` changes that diagnostic grid only. Oversampling an SGF
record does not add temporal resolution. The event curves and crest values
are native samples, not interpolation-grid maxima. No parabolic peak fitting,
amplitude filtering or smoothing is applied.

## Crest-height statistics

Each complete upcrossing wave contributes **one** crest: its largest sampled
elevation minus that record's own time-weighted mean. Report:

* complete-wave count, mean/median crest height, P90/P95/P99/P99.9 and maximum;
* maximum crest-to-trough wave height;
* absolute maximum elevation above still water and its time;
* mean water level, time-weighted surface standard deviation and `Hm0 = 4 sigma`;
* minimum wave period and count of waves shorter than `Tp/2`, exposing short
  waves/ripples rather than silently filtering them out;
* full-record SGF-versus-HOS temporal agreement on absolute, not demeaned,
  elevations, so mean-level errors remain visible.

The mean and variance integrate the piecewise-linear, timestamp-collapsed
history, rather than averaging native samples. This avoids overweighting wave
phases where adaptive timestepping produces more samples. Native extrema are
preserved separately. Small ripples can change upcrossing counts; wave counts
need not match across records, and distribution Q-Q points are **quantile**
comparisons, not a pairing of individual waves.

The crest exceedance figure uses Weibull plotting positions
`rank/(N+1)` in descending crest order. These are empirical sample statistics,
not extrapolated return levels or fitted extreme-value distributions.
With roughly 1200 complete waves, P99.9 and the maximum have particularly
limited statistical support. No confidence bounds or assumption that all
waves are independent are implied.

## Saved results and existing-data example

Default outputs under `plotting_outputs/`:

* `largest_wave_events.png`: six HOS-selected, unaligned temporal comparisons;
* `crest_height_statistics.png`: exceedance curves, distribution Q-Q,
  characteristic crests and a summary table;
* `temporal_comparison.json`: methods, source/clock metadata, record statistics,
  event definitions, matched crest errors and timing diagnostics;
* `crest_statistics.csv`, `largest_wave_events.csv`, `wave_crests.csv`:
  summary statistics, all matched events and every complete extracted wave;
* `temporal_comparison.npz`: the shared comparison grid/elevations and the
  native crest-height/time populations behind the figures.

For the default **30--1776.99 s** window:

| Record | Complete waves | P99 crest [m] | Maximum crest [m] | Temporal RMSE vs HOS [m] |
|---|---:|---:|---:|---:|
| HOS-NWT | 1187 | 0.12394 | 0.19864 | -- |
| SGF, dz=0.0188 m | 1176 | 0.11786 | 0.16414 | 0.003884 |
| SGF, dz=0.0094 m | 1176 | 0.12224 | 0.18757 | 0.000632 |

Crests above are relative to each record's time-weighted mean. At the largest
HOS event, **t=362.464 s**, the matched coarse/fine SGF crests are **0.16108 m**
and **0.18757 m**, errors of **-18.91%** and **-5.57%**. The coarse record's
own global maximum occurs at another event. Timing is close for this event
but not uniformly so: the second-largest selected HOS crest has a coarse-SGF
lag of about **0.118 s**, versus about **0.008 s** in the fine result.

These are comparisons of existing simulated records, not new simulations or
experimental validation. No additional sampling is needed for these figures.
For future exports, preserve full-precision physical timestamps to avoid
collapsing adaptive samples and to improve peak-time interpretation.

## Tests

```sh
python3 -m unittest test_temporal_analysis
```

Tests cover HOS clock recovery and rounding checks, tied timestamps without
lost native crests, complete upcrossing waves, adaptive-sampling moments,
reference-based separated event ranking, retained amplitude/timing/level errors,
invalid input, and persistent CSV/JSON/NPZ/figure outputs. The plotting test
runs when Matplotlib is installed.
