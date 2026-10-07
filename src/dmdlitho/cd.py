"""Dose crossings and the mean top/bottom CD illustrated in Fig. 4."""

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq
from scipy.special import erfcinv

from .exposure import normalized_dose_profile, normalized_peak_dose
from .parameters import LithographyParameters, ThresholdMode, dose_thresholds, positive_finite


@dataclass(frozen=True)
class CDResult:
    """Widths in µm; None denotes a physically absent crossing or CD."""

    cd_um: float | None
    cd_top_um: float | None
    cd_bottom_um: float | None
    fully_cleared: bool


def _positive_crossing(width: float, radius: float, level: float) -> float | None:
    """Find the positive-side position where E/E0 falls to level."""
    peak = normalized_peak_dose(width, radius)
    # At clearing onset the bottom touches the substrate at x=0. Eight ulps
    # accommodate rounding in erf(erfinv(level)) and dose scaling, not a
    # physical clearing tolerance. A genuinely missing crossing returns None.
    rounding = 8 * np.finfo(float).eps * max(peak, level)
    if level > peak + rounding:
        return None
    if abs(level - peak) <= rounding:
        return 0.0

    def residual(x: float) -> float:
        return normalized_dose_profile(x, width, radius) - level

    # Beyond the written edge, E/E0 <= 0.5*erfc((x-l/2)/(sqrt(2)*w)).
    # Its inverse supplies a dose-dependent physical bracket instead of an
    # arbitrary huge coordinate. One radius of margin makes the bound strict.
    tail_offset = max(0.0, float(erfcinv(2.0 * level)))
    upper = 0.5 * width + radius * (np.sqrt(2.0) * tail_offset + 1.0)
    # The isolated Gaussian dose decreases continuously for x>0, so Brent's
    # bracketing method finds the unique edge without sampling a plotting grid.
    return float(brentq(residual, 0.0, upper, xtol=1e-11 * radius))


def compute_cd(write_width_um: float, parameters: LithographyParameters,
               dose_scale: float = 1.0,
               threshold_mode: ThresholdMode = "reconstructed") -> CDResult:
    """Compute Fig. 4 CD at scaled dose; return absent edges explicitly."""
    width = positive_finite(write_width_um, "write_width_um")
    scale = positive_finite(dose_scale, "dose_scale")
    threshold, clearing = dose_thresholds(parameters.gamma, threshold_mode)
    top_level, bottom_level = threshold / scale, clearing / scale
    if top_level <= 0 or bottom_level <= 0:
        raise ValueError("dose_scale gives unrepresentable crossing levels")
    x_top = _positive_crossing(width, parameters.gaussian_radius_um, top_level)
    x_bottom = _positive_crossing(width, parameters.gaussian_radius_um, bottom_level)
    # If the center has not reacted, neither edge exists. If it has reacted
    # but has not cleared, the top exists but the substrate edge and mean CD
    # are undefined. Never substitute a zero bottom for incomplete clearing.
    top = None if x_top is None else 2.0 * x_top
    if x_bottom is None:
        return CDResult(None, top, None, False)
    # Symmetry places edges at ±x: CD_top=2*x_top, CD_bottom=2*x_bottom.
    # The paper's arithmetic mean is therefore CD=x_top+x_bottom.
    return CDResult(x_top + x_bottom, top, 2.0 * x_bottom, True)


def matching_index(design_width_um: float, printed_cd_um: float) -> float:
    """Section 4.1 matching index = design width / printed CD * 100 (%)."""
    design = positive_finite(design_width_um, "design_width_um")
    printed = positive_finite(printed_cd_um, "printed_cd_um")
    return 100.0 * design / printed
