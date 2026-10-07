"""Simplified positive-resist vertical development, not calibrated kinetics."""

import numpy as np

from .parameters import positive_finite


def normalized_removed_depth(normalized_dose, threshold_ratio: float,
                             clearing_ratio: float, gamma: float):
    """Eq. (6): d/T0 = gamma ln(E/Eth), clipped to [0,1].

    All doses are dimensionless ratios to the same E0. Zero means no
    removed resist; one means fully cleared. Remaining thickness is
    T0*(1-d/T0), avoiding the paper's ambiguous t(x) terminology.
    """
    threshold_ratio = positive_finite(threshold_ratio, "threshold_ratio")
    clearing_ratio = positive_finite(clearing_ratio, "clearing_ratio")
    gamma = positive_finite(gamma, "gamma")
    if clearing_ratio <= threshold_ratio:
        raise ValueError("clearing_ratio must exceed threshold_ratio")
    log_ratio = np.log(clearing_ratio) - np.log(threshold_ratio)
    if not np.isclose(gamma * log_ratio, 1.0, rtol=1e-10, atol=1e-12):
        raise ValueError("thresholds and gamma must satisfy Eq. (5)")
    dose = np.asarray(normalized_dose, dtype=float)
    if not np.all(np.isfinite(dose) & (dose >= 0)):
        raise ValueError("normalized_dose must be finite and nonnegative")
    # Only doses above Eth enter the logarithm. This is a physical reaction
    # cutoff and avoids log(0), rather than introducing an artificial dose floor.
    removed = np.zeros_like(dose)
    active = (dose > threshold_ratio) & (dose < clearing_ratio)
    removed[active] = gamma * (np.log(dose[active]) - np.log(threshold_ratio))
    removed[dose >= clearing_ratio] = 1.0
    removed = np.clip(removed, 0.0, 1.0)
    return float(removed) if dose.ndim == 0 else removed
