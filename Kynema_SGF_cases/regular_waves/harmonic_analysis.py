"""Joint harmonic least-squares fits at a common, target-defined frequency."""

import numpy as np

from reflection.reflection_coeff_komega import resample_uniform


def sample_spacing(times):
    """(median, maximum) output interval of increasing, possibly non-uniform times."""
    times = np.asarray(times, dtype=float)
    if times.ndim != 1 or times.size < 2 or not np.isfinite(times).all():
        raise ValueError("time must contain at least two finite values")
    spacing = np.diff(times)
    if spacing.min() <= 0:
        raise ValueError("time must be strictly increasing")
    return float(np.median(spacing)), float(spacing.max())


def design_matrix(times, frequency, count):
    clock = 2 * np.pi * frequency * (times - times[0])
    return np.column_stack(
        [np.ones(len(times))]
        + [component for h in range(1, count + 1)
           for component in (np.cos(h * clock), np.sin(h * clock))]
    )


def estimate_frequency(times, target, frequency_band, count):
    """Refine the target FFT peak by variable-projection least squares.

    The FFT peak (initial guess only) uses probes resampled to uniform time if the
    outputs are non-uniform; the least-squares refinement uses the actual times.
    """
    low, high = frequency_band
    probes = target[:, np.unique(np.linspace(0, target.shape[1] - 1, min(16, target.shape[1])).astype(int))]
    if not np.any(probes - probes.mean(axis=0)):
        raise ValueError("Target has no oscillatory signal for a fundamental-frequency fit")
    grid, uniform_probes, _ = resample_uniform(times, probes)
    dt = float(grid[1] - grid[0])
    centered = uniform_probes - uniform_probes.mean(axis=0)
    frequencies = np.fft.rfftfreq(len(grid), dt)
    power = np.sum(abs(np.fft.rfft(centered * np.hanning(len(grid))[:, None], axis=0))**2, axis=1)
    in_band = np.flatnonzero((frequencies >= low) & (frequencies <= high))
    if not len(in_band):
        raise ValueError("No FFT frequencies in the fundamental band; use a longer record")
    peak = frequencies[in_band[np.argmax(power[in_band])]]
    width = 1 / (len(grid) * dt)

    def objective(frequency):
        matrix = design_matrix(times, frequency, count)
        coefficients = np.linalg.lstsq(matrix, probes, rcond=None)[0]
        return float(np.sum((probes - np.einsum("ij,jk->ik", matrix, coefficients))**2))

    grid = np.linspace(max(low, peak - width), min(high, peak + width), 17)
    best = int(np.argmin([objective(frequency) for frequency in grid]))
    left, right = grid[max(0, best - 1)], grid[min(len(grid) - 1, best + 1)]
    ratio = (np.sqrt(5) - 1) / 2
    c, d = right - ratio * (right - left), left + ratio * (right - left)
    fc, fd = objective(c), objective(d)
    while right - left > 1e-10:
        if fc < fd:
            right, d, fd = d, c, fc
            c = right - ratio * (right - left)
            fc = objective(c)
        else:
            left, c, fc = c, d, fd
            d = left + ratio * (right - left)
            fd = objective(d)
    frequency = float((left + right) / 2)
    if min(frequency - low, high - frequency) < 1e-7:
        raise ValueError("Target frequency fit reached a band edge; widen the fundamental band")
    return frequency


def harmonic_analysis(times, eta, target, frequency_band=(0.4, 0.85), *,
                      count=5, frequency=None, phase_floor_ratio=1e-3):
    """Fit eta = mean + Re(sum(C_h exp(-i h omega (t-t0)))) at every x.

    Positive arg(C_signal*conj(C_target)) is a temporal lag at a fixed x.
    No time/phase alignment or independent signal-frequency fitting is applied.
    """
    times = np.asarray(times, dtype=float)
    typical_dt, max_dt = sample_spacing(times)
    eta, target = np.asarray(eta, dtype=float), np.asarray(target, dtype=float)
    if (eta.ndim != 2 or eta.shape != target.shape or eta.shape[0] != len(times)
            or eta.shape[1] < 1 or not np.isfinite(eta).all() or not np.isfinite(target).all()):
        raise ValueError("Signal and target must be finite matching arrays with shape (time, point)")
    if not isinstance(count, int) or count < 1 or len(times) < 2 * count + 2:
        raise ValueError("Insufficient samples for the requested joint harmonic fit")
    low, high = frequency_band
    if not np.isfinite([low, high]).all() or not 0 < low < high:
        raise ValueError("Require a finite, positive, increasing fundamental-frequency band")
    if not 0 < phase_floor_ratio < 1:
        raise ValueError("phase_floor_ratio must lie strictly between zero and one")
    estimate_from_target = frequency is None
    if estimate_from_target:
        # The FFT peak search runs on data resampled at the median interval.
        if count * high >= 0.5 / typical_dt:
            raise ValueError("Highest fitted harmonic would reach Nyquist within the search band")
        frequency = estimate_frequency(times, target, frequency_band, count)
    # Nyquist uses the median output interval. Under adaptive dt the outputs are
    # irregular (occasional gaps approach twice the interval); the least-squares fit
    # does not alias on irregular times, and an ill-conditioned fit raises below.
    if not np.isfinite(frequency) or not low < frequency < high or count * frequency >= 0.5 / typical_dt:
        raise ValueError("Fundamental must be inside the band and all harmonics below Nyquist")
    matrix = design_matrix(times, frequency, count)
    coefficients, _, rank, singular = np.linalg.lstsq(
        matrix, np.column_stack((eta, target)), rcond=None
    )
    if rank != matrix.shape[1] or singular[0] / singular[-1] > 1e8:
        raise ValueError("Harmonic fit is ill-conditioned; use a longer record or fewer harmonics")
    nx = eta.shape[1]
    signal_c = coefficients[1::2, :nx] + 1j * coefficients[2::2, :nx]
    target_c = coefficients[1::2, nx:] + 1j * coefficients[2::2, nx:]
    signal_a, target_a = abs(signal_c), abs(target_c)
    scale = float(np.sqrt(np.mean(target_a[0]**2)))
    if scale <= np.finfo(float).eps * max(float(np.max(abs(target))), np.finfo(float).tiny) * 100:
        raise ValueError("Target fundamental is too weak to define harmonic amplitudes and phase")
    floor = phase_floor_ratio * scale
    valid = (target_a >= floor) & (signal_a >= floor)
    phase = (np.angle(signal_c * target_c.conj()) + np.pi) % (2 * np.pi) - np.pi
    phase[~valid] = np.nan
    orders = np.arange(1, count + 1)
    reports = []
    for index, h in enumerate(orders):
        signal_rms = float(np.sqrt(np.mean(signal_a[index]**2)))
        target_rms = float(np.sqrt(np.mean(target_a[index]**2)))
        mask = valid[index]
        reference_resolved = target_rms >= floor
        mean_phase = rms_phase = max_phase = concentration = None
        status = "weak_reference_or_signal"
        if mask.any():
            weights = target_a[index, mask]**2
            angles = phase[index, mask]
            vector = np.sum(weights * np.exp(1j * angles)) / weights.sum()
            concentration = float(abs(vector))
            rms_phase = float(np.rad2deg(np.sqrt(np.sum(weights * angles**2) / weights.sum())))
            max_phase = float(np.rad2deg(np.max(abs(angles))))
            status = "resolved" if concentration > 1e-8 else "incoherent_phase"
            if status == "resolved":
                mean_phase = float(np.rad2deg(np.angle(vector)))
        reports.append({
            "order": int(h),
            "frequency_hz": float(h * frequency),
            "signal_amplitude_rms_m": signal_rms,
            "target_amplitude_rms_m": target_rms,
            "signal_amplitude_over_target_fundamental": signal_rms / scale,
            "target_amplitude_over_target_fundamental": target_rms / scale,
            "relative_amplitude_error": signal_rms / target_rms - 1 if reference_resolved else None,
            "phase_status": status,
            "phase_valid_points": int(mask.sum()),
            "phase_error_circular_mean_deg": mean_phase,
            "phase_error_rms_deg": rms_phase,
            "phase_error_max_abs_deg": max_phase,
            "phase_resultant_length": concentration,
            "equivalent_time_lag_s": mean_phase / (360 * h * frequency) if mean_phase is not None else None,
        })
    # Explicit contraction avoids false floating-point flags in some BLAS backends.
    fitted = np.einsum("ij,jk->ik", matrix, coefficients)
    signal_residual = np.sqrt(np.mean((eta - fitted[:, :nx])**2, axis=0))
    target_residual = np.sqrt(np.mean((target - fitted[:, nx:])**2, axis=0))
    summary = {
        "fundamental_frequency_hz": frequency,
        "frequency_source": "target-only joint harmonic fit" if estimate_from_target else "specified common frequency",
        "fit_duration_s": float(times[-1] - times[0] + typical_dt),
        "periods_in_record": float((times[-1] - times[0] + typical_dt) * frequency),
        "short_record": bool((times[-1] - times[0] + typical_dt) * frequency < 10),
        "median_sample_interval_s": typical_dt,
        "max_sample_interval_s": max_dt,
        "phase_convention": "wrapped signal minus target; positive means temporal lag",
        "phase_floor_ratio": phase_floor_ratio,
        "phase_amplitude_floor_m": floor,
        "time_origin_s": float(times[0]),
        "signal_fit_residual_rms_m": float(np.sqrt(np.mean(signal_residual**2))),
        "target_fit_residual_rms_m": float(np.sqrt(np.mean(target_residual**2))),
        "harmonics": reports,
    }
    return summary, {
        "harmonic_orders": orders,
        "harmonic_frequencies_hz": orders * frequency,
        "harmonic_signal_coefficient_m": signal_c,
        "harmonic_target_coefficient_m": target_c,
        "harmonic_signal_amplitude_m": signal_a,
        "harmonic_target_amplitude_m": target_a,
        "harmonic_phase_error_rad": phase,
        "harmonic_phase_valid": valid,
        "harmonic_signal_residual_rms_m": signal_residual,
        "harmonic_target_residual_rms_m": target_residual,
    }


def plot_harmonics(summary, arrays, output):
    """Plot case harmonic amplitudes and local, wrapped phase errors."""
    import matplotlib.pyplot as plt

    harmonics = summary["harmonic_analysis"]["harmonics"]
    orders = arrays["harmonic_orders"]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4), layout="constrained")
    axes[0].semilogy(orders, [row["signal_amplitude_rms_m"] for row in harmonics], "o-", label="Computed")
    axes[0].semilogy(orders, [row["target_amplitude_rms_m"] for row in harmonics], "s--", label="Target")
    axes[0].set(xlabel="Harmonic order", ylabel="Spatial RMS of amplitude [m]", xticks=orders)
    axes[0].legend()
    axes[1].plot(orders, [row["phase_error_circular_mean_deg"] for row in harmonics], "o-")
    axes[1].set(xlabel="Harmonic order", ylabel="Circular-mean phase error [deg]", xticks=orders)
    for index, h in enumerate(orders):
        mask = arrays["harmonic_phase_valid"][index]
        axes[2].plot(arrays["x_m"][mask], np.rad2deg(arrays["harmonic_phase_error_rad"][index, mask]),
                     label=f"h={h}")
    axes[2].set(xlabel="x [m]", ylabel="Local wrapped phase error [deg]")
    axes[2].legend(ncol=2)
    for ax in axes:
        ax.grid(True, alpha=0.25)
    axes[1].axhline(0, color="gray", linewidth=0.7)
    axes[2].axhline(0, color="gray", linewidth=0.7)
    fig.suptitle(f"{summary['case']}: harmonics and phase (positive = lag)")
    fig.savefig(output, dpi=180)
    plt.close(fig)
