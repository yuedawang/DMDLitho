import numpy as np
import pytest
from scipy.integrate import quad

from dmdlitho import normalized_dose_profile
from dmdlitho.scanning import (DMDScanParameters, build_address_cells,
                               discrete_normalized_dose_profile, profile_error_metrics)


def test_fig3_alignment_geometry():
    scan = DMDScanParameters(writing_period=4)
    assert scan.tilt_angle_rad == pytest.approx(np.arctan(1/4))
    assert scan.address_pitch_um == pytest.approx(25*np.cos(np.arctan(1/4))/4)
    assert scan.redundancy == 192


def test_paper_experimental_scan_condition():
    scan = DMDScanParameters()
    assert scan.tilt_angle_deg == pytest.approx(3.58, abs=0.01)
    assert scan.address_pitch_um == pytest.approx(1.56, abs=0.003)
    assert scan.moving_pitch_um == pytest.approx(24000/15400)
    assert scan.moving_pitch_um == pytest.approx(1.56, abs=0.003)
    assert scan.redundancy == 48
    assert scan.repeated_writing_time_s == pytest.approx(48/15400)
    assert scan.relative_pitch_mismatch < 0.005


def test_units_and_exact_synchronization():
    # Choose f from the geometric pitch, rather than rounding either pitch.
    geometry = DMDScanParameters()
    scan = DMDScanParameters(switching_frequency_hz=24000/geometry.address_pitch_um)
    assert scan.moving_pitch_um == pytest.approx(scan.address_pitch_um)
    assert scan.relative_pitch_mismatch == pytest.approx(0, abs=1e-15)


@pytest.mark.parametrize("name", ["writing_period", "array_columns"])
@pytest.mark.parametrize("bad", [0, -1, 4.5, 16.0, True])
def test_invalid_integer_geometry(name, bad):
    with pytest.raises(ValueError, match="positive integer"):
        DMDScanParameters(**{name: bad})


def test_noninteger_redundancy_is_not_truncated():
    with pytest.raises(ValueError, match="divisible"):
        DMDScanParameters(writing_period=7)


@pytest.mark.parametrize("name", ["beam_pitch_um", "scan_velocity_mm_s", "switching_frequency_hz"])
@pytest.mark.parametrize("bad", [0, -1, np.nan, np.inf])
def test_invalid_physical_scan_inputs(name, bad):
    with pytest.raises(ValueError):
        DMDScanParameters(**{name: bad})


@pytest.mark.parametrize("width,pitch", [(10, 2), (10, 3), (10, 20), (1.1, 0.2)])
def test_centered_cells_preserve_width_and_symmetry(width, pitch):
    centers, weights = build_address_cells(width, pitch)
    np.testing.assert_allclose(centers, -centers[::-1], atol=1e-14)
    np.testing.assert_allclose(weights, weights[::-1], atol=1e-14)
    assert np.all(weights > 0) and np.all(weights <= pitch)
    assert np.all(np.abs(centers) < width/2)
    assert np.sum(weights) == pytest.approx(width, abs=1e-13)


def test_partial_boundaries_have_actual_midpoints_and_weights():
    # Three full 3 µm cells span [-4.5,4.5]; two 0.5 µm cells close the
    # exact [-5,5] interval, without shifting its width or breaking symmetry.
    centers, weights = build_address_cells(10, 3)
    np.testing.assert_allclose(centers, [-4.75, -3, 0, 3, 4.75])
    np.testing.assert_allclose(weights, [0.5, 3, 3, 3, 0.5])


def test_symmetric_nonnegative_profile_and_scalar_input():
    positions = np.linspace(-17, 17, 801)
    discrete = discrete_normalized_dose_profile(positions, 10, 2, 1.56)
    np.testing.assert_allclose(discrete, discrete_normalized_dose_profile(-positions, 10, 2, 1.56),
                               atol=1e-14, rtol=1e-14)
    assert np.all(discrete >= 0)
    assert discrete_normalized_dose_profile(positions[201], 10, 2, 1.56) == pytest.approx(discrete[201])


def test_single_cell_has_quadrature_amplitude_and_total_integral():
    # A coarse one-cell interval contributes l*g(x). The integral must be l,
    # not 1; this independently detects hidden peak/area normalization.
    positions = np.array([0, 1, 3])
    expected = 1.2 * np.exp(-positions**2/8) / (2*np.sqrt(2*np.pi))
    np.testing.assert_allclose(discrete_normalized_dose_profile(positions, 1.2, 2, 3), expected)
    area = quad(lambda x: discrete_normalized_dose_profile(x, 1.2, 2, 3), -np.inf, np.inf)[0]
    assert area == pytest.approx(1.2, abs=1e-10)


def test_geometry_sequence_reduces_both_profile_errors():
    positions = np.linspace(-17, 17, 6001)
    analytical = normalized_dose_profile(positions, 10, 2)
    errors = []
    for period in [4, 8, 16, 32, 64]:
        pitch = DMDScanParameters(writing_period=period).address_pitch_um
        discrete = discrete_normalized_dose_profile(positions, 10, 2, pitch)
        metrics = profile_error_metrics(discrete, analytical)
        errors.append([metrics.relative_l2_error, metrics.max_absolute_error])
    assert np.all(np.diff(errors, axis=0) < 0)
    assert errors[-1][0] < 0.002 and errors[-1][1] < 0.002


def test_fine_pitch_agreement_and_evaluation_grid_resolution():
    metrics_by_grid = []
    for count in [3001, 6001]:
        positions = np.linspace(-17, 17, count)
        reference = normalized_dose_profile(positions, 10, 2)
        discrete = discrete_normalized_dose_profile(positions, 10, 2, 0.05)
        np.testing.assert_allclose(discrete, reference, atol=2e-5, rtol=0)
        metrics = profile_error_metrics(discrete, reference)
        metrics_by_grid.append(metrics.relative_l2_error)
    assert metrics_by_grid[0] == pytest.approx(metrics_by_grid[1], rel=1e-5)


def test_profile_error_definitions():
    metrics = profile_error_metrics([3, 5], [3, 4])
    assert metrics.relative_l2_error == pytest.approx(1/5)
    assert metrics.max_absolute_error == 1


@pytest.mark.parametrize("discrete,reference", [([], []), ([1], [1, 2]),
                                                 ([[1]], [[1]]), ([np.nan], [1]),
                                                 ([1], [np.inf]), ([1], [0])])
def test_invalid_profile_metric_inputs(discrete, reference):
    with pytest.raises(ValueError):
        profile_error_metrics(discrete, reference)


@pytest.mark.parametrize("bad", [0, -1, np.nan, np.inf])
def test_invalid_discrete_lengths(bad):
    with pytest.raises(ValueError):
        build_address_cells(bad, 1)
    with pytest.raises(ValueError):
        build_address_cells(10, bad)
    with pytest.raises(ValueError):
        discrete_normalized_dose_profile(0, 10, bad, 1)


def test_invalid_discrete_positions():
    with pytest.raises(ValueError):
        discrete_normalized_dose_profile([0, np.nan], 10, 2, 1)
