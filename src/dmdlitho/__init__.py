"""DMDLitho v0.1: isolated-line numerical reproduction of Hur & Seo (2012)."""

from .parameters import LithographyParameters, literal_thresholds, reconstructed_thresholds
from .exposure import normalized_dose_profile, normalized_peak_dose
from .resist import normalized_removed_depth
from .cd import CDResult, compute_cd, matching_index
from .correction import (UnachievableTargetError, clearing_write_width,
                         solve_width_correction, solve_dose_correction, resolution_limits)

__all__ = ["LithographyParameters", "literal_thresholds", "reconstructed_thresholds",
           "normalized_dose_profile", "normalized_peak_dose", "normalized_removed_depth",
           "CDResult", "compute_cd", "matching_index", "UnachievableTargetError",
           "clearing_write_width", "solve_width_correction", "solve_dose_correction",
           "resolution_limits"]
