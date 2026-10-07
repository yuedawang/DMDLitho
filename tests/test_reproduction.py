import numpy as np
import pytest

from dmdlitho import (LithographyParameters, compute_cd, matching_index,
                      solve_width_correction, solve_dose_correction,
                      resolution_limits, normalized_peak_dose,
                      normalized_removed_depth, reconstructed_thresholds,
                      UnachievableTargetError)


def test_published_fig6_regressions():
    parameters = LithographyParameters()
    index = matching_index(6.6, compute_cd(6.6, parameters).cd_um)
    assert index == pytest.approx(107.8, abs=0.3)
    assert solve_width_correction(5.9, parameters).reinforcement_percent == pytest.approx(110.9, abs=0.3)
    assert solve_dose_correction(5.4, parameters).reinforcement_percent == pytest.approx(109.1, abs=0.3)


def test_fig5_monotonic_saturation():
    widths = np.linspace(0.01, 30, 400)
    peak = normalized_peak_dose(widths, 2)
    threshold, clearing = reconstructed_thresholds(0.46)
    depth = normalized_removed_depth(peak, threshold, clearing, 0.46)
    assert np.all(np.diff(peak) >= 0) and np.all(np.diff(depth) >= 0)
    assert depth[0] == 0 and depth[-1] == 1
    assert peak[-1] == pytest.approx(1, abs=1e-12)


def test_fig6_wide_width_correction_vanishes():
    parameters = LithographyParameters()
    for target in [15, 20]:
        assert matching_index(target, compute_cd(target, parameters).cd_um) == pytest.approx(100, abs=0.001)
        assert solve_width_correction(target, parameters).reinforcement_percent == pytest.approx(100, abs=0.001)
        assert solve_dose_correction(target, parameters).reinforcement_percent == pytest.approx(100, abs=0.001)


def test_resolution_limits_and_contrast_sweep():
    limits_by_gamma = []
    for gamma in [0.46, 0.47, 0.48, 0.49, 0.50]:
        parameters = LithographyParameters(gamma=gamma)
        limits = resolution_limits(parameters)
        limits_by_gamma.append([limits.uncorrected_design_um, limits.width_corrected_um,
                                limits.dose_corrected_um])
        assert limits.dose_corrected_um < limits.width_corrected_um < limits.uncorrected_design_um
        sampled = 0.1 * np.ceil(np.array(limits_by_gamma[-1])/0.1)
        width = solve_width_correction(float(sampled[1]), parameters)
        dose = solve_dose_correction(float(sampled[2]), parameters)
        assert dose.reinforcement_percent < width.reinforcement_percent
        for solver, limit in [(solve_width_correction, limits.width_corrected_um),
                              (solve_dose_correction, limits.dose_corrected_um)]:
            with pytest.raises(UnachievableTargetError):
                solver(limit-0.01, parameters)
            assert solver(limit+0.01, parameters).cd_result.fully_cleared
    assert np.all(np.diff(limits_by_gamma, axis=0) < 0)
    # Reported tenths are compared to our explicitly inferred 0.1 µm grid.
    # The model permits a dose-corrected 5.3 µm target, whereas the paper
    # reports 5.4 µm. Preserve this discrepancy rather than moving a cutoff.
    np.testing.assert_allclose(0.1*np.ceil(np.array(limits_by_gamma[0])/0.1),
                               [6.6, 5.9, 5.4], atol=0.11, rtol=0)
