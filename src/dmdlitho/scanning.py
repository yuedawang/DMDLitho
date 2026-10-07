"""Paper-aligned scan geometry and a 1D midpoint approximation to Eq. (4).

This module calculates geometry and normalized spatial accumulation. It does
not simulate DMD frames, a two-dimensional array, or a physical stage.
"""

from dataclasses import dataclass
from numbers import Integral

import numpy as np

from .parameters import positive_finite


def _positive_integer(value: int, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, Integral) or value <= 0:
        raise ValueError(f"{name} must be a positive integer")


@dataclass(frozen=True)
class DMDScanParameters:
    """Projected-array alignment parameters; p in µm, v in mm/s, f in Hz.

    Defaults reconstruct the Section 4.2 scan condition, not its corner
    correction. N must be divisible by A so R counts whole repeated writes.
    """

    beam_pitch_um: float = 25.0
    writing_period: int = 16
    array_columns: int = 768
    scan_velocity_mm_s: float = 24.0
    switching_frequency_hz: float = 15400.0

    def __post_init__(self) -> None:
        for name in ("beam_pitch_um", "scan_velocity_mm_s", "switching_frequency_hz"):
            positive_finite(getattr(self, name), name)
        _positive_integer(self.writing_period, "writing_period")
        _positive_integer(self.array_columns, "array_columns")
        if self.array_columns % self.writing_period != 0:
            raise ValueError("array_columns must be divisible by writing_period for integer redundancy")

    @property
    def tilt_angle_rad(self) -> float:
        """theta=atan(1/A), radians, for the Fig. 3 row/column alignment."""
        # Aligning (i,j) and (i+1,j+A) with the scan direction gives transverse
        # and longitudinal offsets in the ratio 1:A. This particular square-
        # lattice construction gives tan(theta)=1/A, not a universal DMD tilt.
        # A=1 is its formal 45-degree limit; the paper's R discussion uses
        # theta<45 degrees (A>1), as do all reproduction cases in our scripts.
        return float(np.arctan(1.0 / self.writing_period))

    @property
    def tilt_angle_deg(self) -> float:
        """Array tilt in degrees for readable geometry reporting."""
        return float(np.degrees(self.tilt_angle_rad))

    @property
    def address_pitch_um(self) -> float:
        """Eq. (2): transverse separation of adjacent address lines, µm."""
        # Tilting interleaves A address lines within the projected row pitch.
        # The cosine projects p normal to the scan before division by A.
        projected_row_pitch_um = self.beam_pitch_um * np.cos(self.tilt_angle_rad)
        return float(projected_row_pitch_um / self.writing_period)

    @property
    def moving_pitch_um(self) -> float:
        """Stage displacement per DMD switching interval, v/f in µm."""
        scan_velocity_um_s = 1000.0 * self.scan_velocity_mm_s
        return scan_velocity_um_s / self.switching_frequency_hz

    @property
    def relative_pitch_mismatch(self) -> float:
        """Absolute synchronization mismatch |p_mov-p_add|/p_add."""
        # p_mov is along the scan, p_add is across it. Approximate equality
        # makes the address spacing comparable in the two directions; it
        # does not mean our one-dimensional model simulates stage motion.
        return abs(self.moving_pitch_um - self.address_pitch_um) / self.address_pitch_um

    @property
    def redundancy(self) -> int:
        """R=N/A repeated writes in the adopted alignment geometry."""
        # Divisibility is checked in __post_init__; never truncate a fractional
        # write count. Changing A at fixed N also changes physical exposure.
        return int(self.array_columns // self.writing_period)

    @property
    def repeated_writing_time_s(self) -> float:
        """Eq. (3) exposure-time factor R/f (seconds), not normalized dose."""
        # Multiplying intensity I(x) by R/f gives the physical single-address
        # accumulated dose e(x). Without an intensity calibration this time
        # factor alone is not a dose in energy per area.
        return self.redundancy / self.switching_frequency_hz


def build_address_cells(width_um: float, address_pitch_um: float) -> tuple[np.ndarray, np.ndarray]:
    """Return symmetric midpoint centers and cell weights (both in µm).

    Full cells have width p_add. Any remainder is split equally between
    the two feature boundaries; those partial cells use their actual
    midpoint and width. A feature narrower than one pitch uses one cell.
    """
    width = positive_finite(width_um, "width_um")
    pitch = positive_finite(address_pitch_um, "address_pitch_um")
    cell_count_ratio = width / pitch
    if not np.isfinite(cell_count_ratio) or cell_count_ratio > np.iinfo(np.intp).max:
        raise ValueError("address_pitch_um cannot represent the requested cell count")
    full_count = int(np.floor(cell_count_ratio))
    if full_count == 0:
        return np.array([0.0]), np.array([width])

    # Place full cells symmetrically about x=0, retaining the requested pitch.
    # Their centers are separated by pitch and their union spans full_span.
    full_span = full_count * pitch
    centers = (np.arange(full_count, dtype=float) + 0.5 - 0.5 * full_count) * pitch
    weights = np.full(full_count, pitch)
    remainder = width - full_span
    if remainder > 0:
        # Cutting only the rightmost cell would break feature symmetry. Two
        # equal boundary cells preserve [-l/2,l/2] without snapping l to a
        # lattice. Their centers need not be spaced by pitch from a neighbor:
        # this is a feature-centered quadrature convention, not hardware-grid
        # quantization or a fixed global DMD address phase.
        boundary_center = 0.25 * (width + full_span)
        centers = np.concatenate(([-boundary_center], centers, [boundary_center]))
        weights = np.concatenate(([0.5 * remainder], weights, [0.5 * remainder]))
    return centers, weights


def discrete_normalized_dose_profile(x_um, width_um: float, gaussian_radius_um: float,
                                     address_pitch_um: float):
    """Dimensionless weighted Gaussian sum approaching Eq. (4); lengths in µm."""
    radius = positive_finite(gaussian_radius_um, "gaussian_radius_um")
    centers, weights = build_address_cells(width_um, address_pitch_um)
    positions = np.asarray(x_um, dtype=float)
    if not np.all(np.isfinite(positions)):
        raise ValueError("x_um must contain finite positions")
    dose = np.zeros_like(positions)
    gaussian_normalization_um = np.sqrt(2.0 * np.pi) * radius

    # Eq. (1) uses the same exp(-x^2/(2*w^2)) convention as v0.1. The kernel
    # g has unit integral and units 1/µm. Eq. (4)/E0 is the integral of g over
    # the unit-amplitude written interval. Each cell therefore contributes
    # weight_um*g(x-center): weight_um is pitch for a full cell, and its
    # actual width for a partial cell. Omitting this weight changes amplitude
    # with address density and does not approximate the convolution integral.
    for center_um, weight_um in zip(centers, weights):
        distance_in_radii = (positions - center_um) / radius
        gaussian_kernel = np.exp(-0.5 * distance_in_radii**2) / gaussian_normalization_um
        dose += weight_um * gaussian_kernel

    # Eq. (3)'s physical R/f scale is absorbed by E0 in the normalized
    # comparison. It is intentionally not multiplied into this spatial sum.
    # Changing A changes both p_add and real repeated-writing dose; the study
    # factors out the overall exposure scale to isolate spatial discretization.
    # As maximum cell width decreases, this midpoint sum tends to Eq. (4).
    # No peak normalization, clipping, or fitted amplitude correction is used:
    # coarse grids can produce peaks above the continuous saturation reference.
    return float(dose) if positions.ndim == 0 else dose


@dataclass(frozen=True)
class ProfileErrors:
    """Dimensionless numerical profile discrepancies, not CD errors or EPE."""

    relative_l2_error: float
    max_absolute_error: float


def profile_error_metrics(discrete_profile, analytical_profile) -> ProfileErrors:
    """Compare finite 1D profiles evaluated on the same uniform x grid."""
    discrete = np.asarray(discrete_profile, dtype=float)
    analytical = np.asarray(analytical_profile, dtype=float)
    if discrete.ndim != 1 or discrete.size == 0 or discrete.shape != analytical.shape:
        raise ValueError("profiles must be nonempty one-dimensional arrays of equal shape")
    if not np.all(np.isfinite(discrete)) or not np.all(np.isfinite(analytical)):
        raise ValueError("profiles must contain finite values")
    reference_norm = float(np.linalg.norm(analytical))
    if reference_norm == 0:
        raise ValueError("analytical_profile must have a nonzero L2 norm")
    difference = discrete - analytical
    # Uniform evaluation spacing cancels from the ratio of L2 norms. These
    # errors compare dose shapes on a finite observation window; they measure
    # neither resist edge displacement (EPE) nor printed critical dimension.
    return ProfileErrors(float(np.linalg.norm(difference) / reference_norm),
                         float(np.max(np.abs(difference))))
