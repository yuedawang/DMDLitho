"""Reported physical parameters and explicitly inferred dose normalization."""

from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy.special import expit

ThresholdMode = Literal["literal", "reconstructed"]


def positive_finite(value: float, name: str) -> float:
    """Validate a positive scalar physical input."""
    value = float(value)
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return value


@dataclass(frozen=True)
class LithographyParameters:
    """Paper defaults: contrast gamma, beam radius w and thickness T0 (µm)."""

    gamma: float = 0.46
    gaussian_radius_um: float = 2.0
    resist_thickness_um: float = 1.0

    def __post_init__(self) -> None:
        for name in ("gamma", "gaussian_radius_um", "resist_thickness_um"):
            positive_finite(getattr(self, name), name)


def literal_thresholds(gamma: float) -> tuple[float, float]:
    """Return dimensionless Eth/E0=0.1 and Ecl/E0 using Eq. (5)."""
    gamma = positive_finite(gamma, "gamma")
    # 0.1 is the rounded value reported in Section 4.1, not an inferred fit.
    log_clearing = np.log(0.1) + 1.0 / gamma
    if log_clearing >= np.log(np.finfo(float).max):
        raise ValueError("gamma gives an unrepresentable literal clearing dose")
    if log_clearing <= np.log(0.1):
        raise ValueError("gamma cannot resolve distinct thresholds in float64")
    clearing = float(np.exp(log_clearing))
    if clearing <= 0.1:
        raise ValueError("gamma cannot resolve distinct thresholds in float64")
    return 0.1, clearing


def reconstructed_thresholds(gamma: float) -> tuple[float, float]:
    """Return inferred Eth/E0 and Ecl/E0 satisfying Eq. (5) and zero bias."""
    gamma = positive_finite(gamma, "gamma")
    # REPRODUCTION INFERENCE, NOT an equation stated by Hur & Seo:
    # Eth/E0 + Ecl/E0 = 1 makes the two wide-line edge offsets cancel.
    # Combining it with the reported Eq. (5), gamma=1/ln(Ecl/Eth),
    # gives Eth/E0=1/(1+exp(1/gamma)). The reported gamma is 0.46;
    # the resulting Eth/E0 (~0.1021) is close to the reported rounded 0.1.
    # expit evaluates the logistic expression without overflowing exp.
    threshold = float(expit(-1.0 / gamma))
    clearing = 1.0 - threshold
    if threshold <= 0 or clearing >= 1 or threshold >= clearing:
        raise ValueError("gamma cannot resolve distinct thresholds in float64")
    return threshold, clearing


def dose_thresholds(gamma: float, mode: ThresholdMode) -> tuple[float, float]:
    """Select a threshold convention explicitly; never mix normalizations."""
    if mode == "literal":
        return literal_thresholds(gamma)
    if mode == "reconstructed":
        return reconstructed_thresholds(gamma)
    raise ValueError("threshold_mode must be 'literal' or 'reconstructed'")
