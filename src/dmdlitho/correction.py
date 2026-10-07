"""Two physically distinct scalar OPC inverse problems, with no fitting."""

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq
from scipy.special import erfinv

from .cd import CDResult, compute_cd
from .exposure import normalized_peak_dose
from .parameters import LithographyParameters, ThresholdMode, dose_thresholds, positive_finite


class UnachievableTargetError(ValueError):
    """The requested target is below the fully-cleared model's limit."""


@dataclass(frozen=True)
class WidthCorrection:
    write_width_um: float
    reinforcement_percent: float
    cd_result: CDResult


@dataclass(frozen=True)
class DoseCorrection:
    dose_scale: float
    reinforcement_percent: float
    cd_result: CDResult


def clearing_write_width(parameters: LithographyParameters,
                         threshold_mode: ThresholdMode = "reconstructed") -> float:
    """Eq. (7) inverse: minimum written width (µm) clearing at unit dose."""
    _, clearing = dose_thresholds(parameters.gamma, threshold_mode)
    if clearing >= 1:
        raise UnachievableTargetError("unit saturation dose cannot reach Ecl")
    return float(2 * np.sqrt(2.0) * parameters.gaussian_radius_um * erfinv(clearing))


def solve_width_correction(target_width_um: float, parameters: LithographyParameters,
                           threshold_mode: ThresholdMode = "reconstructed") -> WidthCorrection:
    """Change written geometry at fixed unit dose to print the target (µm)."""
    target = positive_finite(target_width_um, "target_width_um")
    lower = clearing_write_width(parameters, threshold_mode)

    def residual(width: float) -> float:
        return compute_cd(width, parameters, threshold_mode=threshold_mode).cd_um - target

    # At this lower bound the bottom touches x=0. Smaller written widths do
    # not clear, so their absent CD cannot be fed into an inverse solver.
    if residual(lower) > 1e-10 * parameters.gaussian_radius_um:
        raise UnachievableTargetError("target is below the width-correction clearing limit")
    if abs(residual(lower)) <= 1e-10 * parameters.gaussian_radius_um:
        width = lower
    else:
        upper = max(target, lower) + 2 * parameters.gaussian_radius_um
        # Increasing the written interval increases dose at every position,
        # hence both edges move outward. Expand by one beam diameter if needed.
        for _ in range(64):
            if residual(upper) >= 0:
                break
            upper += 2 * parameters.gaussian_radius_um
        else:
            raise ValueError("could not bracket width correction at float64 precision")
        width = float(brentq(residual, lower, upper, xtol=1e-11 * parameters.gaussian_radius_um))
    result = compute_cd(width, parameters, threshold_mode=threshold_mode)
    return WidthCorrection(width, 100.0 * width / target, result)


def solve_dose_correction(target_width_um: float, parameters: LithographyParameters,
                          threshold_mode: ThresholdMode = "reconstructed") -> DoseCorrection:
    """Keep written width equal to target; scale dose to print that width."""
    target = positive_finite(target_width_um, "target_width_um")
    _, clearing = dose_thresholds(parameters.gamma, threshold_mode)
    peak = normalized_peak_dose(target, parameters.gaussian_radius_um)
    lower = clearing / peak

    def residual(scale: float) -> float:
        return compute_cd(target, parameters, scale, threshold_mode).cd_um - target

    # This dose is the first to reach the substrate at the center. If even
    # that profile's mean CD exceeds the target, increasing dose only widens
    # it further: a fully-cleared solution is physically impossible.
    if residual(lower) > 1e-10 * parameters.gaussian_radius_um:
        raise UnachievableTargetError("target is below the dose-correction clearing limit")
    if abs(residual(lower)) <= 1e-10 * parameters.gaussian_radius_um:
        scale = lower
    else:
        upper = max(1.0, 2.0 * lower)
        # Dose is unbounded in this simplified model; doubling provides a
        # transparent bracket for the monotonically outward-moving edges.
        for _ in range(64):
            if residual(upper) >= 0:
                break
            upper *= 2.0
        else:
            raise ValueError("could not bracket dose correction at float64 precision")
        scale = float(brentq(residual, lower, upper, xtol=1e-12))
    return DoseCorrection(scale, 100.0 * scale,
                          compute_cd(target, parameters, scale, threshold_mode))


@dataclass(frozen=True)
class ResolutionLimits:
    """Continuous limits in µm; uncorrected design and printed CD differ."""

    uncorrected_design_um: float
    uncorrected_printed_cd_um: float
    width_corrected_um: float
    dose_corrected_um: float


def resolution_limits(parameters: LithographyParameters,
                      threshold_mode: ThresholdMode = "reconstructed") -> ResolutionLimits:
    """Find clearing-onset limits without assuming a width sampling pitch."""
    onset_width = clearing_write_width(parameters, threshold_mode)
    onset_cd = compute_cd(onset_width, parameters, threshold_mode=threshold_mode).cd_um
    threshold, clearing = dose_thresholds(parameters.gamma, threshold_mode)

    def dose_onset_residual(width: float) -> float:
        # At minimum clearing dose the bottom edge is zero. The dose-limited
        # target is where the remaining top half-width equals the design width.
        # For very narrow lines the limiting Gaussian top half-width is
        # w*sqrt(2*ln(Ecl/Eth)), furnishing a positive physical lower bound.
        scale = clearing / normalized_peak_dose(width, parameters.gaussian_radius_um)
        return compute_cd(width, parameters, scale, threshold_mode).cd_um - width

    lower = 0.5 * parameters.gaussian_radius_um * np.sqrt(2 * np.log(clearing / threshold))
    upper = onset_width
    while dose_onset_residual(upper) > 0:
        upper *= 2
    dose_limit = float(brentq(dose_onset_residual, lower, upper,
                             xtol=1e-11 * parameters.gaussian_radius_um))
    return ResolutionLimits(onset_width, onset_cd, onset_cd, dose_limit)
