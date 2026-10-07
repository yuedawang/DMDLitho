"""Gaussian accumulated exposure for an isolated one-dimensional line."""

import numpy as np
from scipy.special import erf, erfc

from .parameters import positive_finite


def normalized_dose_profile(x_um, width_um: float, gaussian_radius_um: float):
    """Eq. (4): E(x)/E0, dimensionless; all input lengths are in µm."""
    width_um = positive_finite(width_um, "width_um")
    gaussian_radius_um = positive_finite(gaussian_radius_um, "gaussian_radius_um")
    x = np.asarray(x_um, dtype=float)
    if not np.all(np.isfinite(x)):
        raise ValueError("x_um must contain finite positions")
    # Integrating the projected Gaussian over the written interval produces
    # an erf difference. Reflection makes the isolated-line symmetry explicit.
    sqrt_two_w = np.sqrt(2.0) * gaussian_radius_um
    left_argument = (np.abs(x) + 0.5 * width_um) / sqrt_two_w
    right_argument = (np.abs(x) - 0.5 * width_um) / sqrt_two_w
    normalized_dose = 0.5 * (erf(left_argument) - erf(right_argument))
    # Outside the line, two erfs approach 1. The equivalent erfc difference
    # preserves small tail doses needed to bracket physical crossings.
    normalized_dose = np.where(
        right_argument > 0,
        0.5 * (erfc(right_argument) - erfc(left_argument)),
        normalized_dose,
    )
    return float(normalized_dose) if x.ndim == 0 else normalized_dose


def normalized_peak_dose(width_um, gaussian_radius_um: float):
    """Eq. (7): E(0)/E0 = erf(l/(2 sqrt(2) w)); lengths in µm."""
    gaussian_radius_um = positive_finite(gaussian_radius_um, "gaussian_radius_um")
    widths = np.asarray(width_um, dtype=float)
    if not np.all(np.isfinite(widths) & (widths > 0)):
        raise ValueError("width_um must contain finite positive widths")
    peak = erf(widths / (2.0 * np.sqrt(2.0) * gaussian_radius_um))
    return float(peak) if widths.ndim == 0 else peak
