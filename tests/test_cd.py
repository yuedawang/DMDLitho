import numpy as np
import pytest

from dmdlitho import (LithographyParameters, compute_cd, matching_index,
                      normalized_dose_profile, reconstructed_thresholds,
                      clearing_write_width, solve_width_correction,
                      solve_dose_correction, UnachievableTargetError)


def test_edges_satisfy_physical_dose_crossings():
    parameters = LithographyParameters()
    threshold, clearing = reconstructed_thresholds(parameters.gamma)
    result = compute_cd(8, parameters, 1.1)
    assert result.fully_cleared
    assert 0 < result.cd_bottom_um < result.cd_um < result.cd_top_um
    assert result.cd_um == pytest.approx((result.cd_top_um + result.cd_bottom_um)/2)
    assert 1.1 * normalized_dose_profile(result.cd_top_um/2, 8, 2) == pytest.approx(threshold)
    assert 1.1 * normalized_dose_profile(result.cd_bottom_um/2, 8, 2) == pytest.approx(clearing)


def test_no_reaction_incomplete_clearing_and_onset():
    parameters = LithographyParameters()
    no_reaction = compute_cd(0.1, parameters)
    assert no_reaction == type(no_reaction)(None, None, None, False)
    incomplete = compute_cd(4, parameters)
    assert incomplete.cd_top_um > 0
    assert incomplete.cd_um is None and incomplete.cd_bottom_um is None
    assert not incomplete.fully_cleared
    onset = compute_cd(clearing_write_width(parameters), parameters)
    assert onset.fully_cleared and onset.cd_bottom_um == 0
    assert onset.cd_um == pytest.approx(onset.cd_top_um/2)


def test_wide_feature_zero_bias_and_literal_difference():
    parameters = LithographyParameters()
    assert compute_cd(30, parameters).cd_um == pytest.approx(30, abs=1e-8)
    assert abs(compute_cd(30, parameters, threshold_mode="literal").cd_um - 30) > 0.01


@pytest.mark.parametrize("bad", [0, -1, np.nan, np.inf])
def test_invalid_inputs(bad):
    parameters = LithographyParameters()
    with pytest.raises(ValueError):
        compute_cd(bad, parameters)
    with pytest.raises(ValueError):
        compute_cd(8, parameters, bad)
    with pytest.raises(ValueError):
        matching_index(8, bad)
    with pytest.raises(ValueError):
        matching_index(bad, 8)
    for solver in [solve_width_correction, solve_dose_correction]:
        with pytest.raises(ValueError):
            solver(bad, parameters)


def test_modes_and_unachievable_targets():
    parameters = LithographyParameters()
    with pytest.raises(ValueError, match="threshold_mode"):
        compute_cd(8, parameters, threshold_mode="unknown")
    for solver in [solve_width_correction, solve_dose_correction]:
        with pytest.raises(UnachievableTargetError):
            solver(1, parameters)
    with pytest.raises(UnachievableTargetError):
        clearing_write_width(LithographyParameters(gamma=0.3), "literal")


@pytest.mark.parametrize("mode", ["literal", "reconstructed"])
@pytest.mark.parametrize("target", [6, 7, 12, 20])
def test_distinct_correction_problems_reach_target(target, mode):
    parameters = LithographyParameters()
    width = solve_width_correction(target, parameters, mode)
    dose = solve_dose_correction(target, parameters, mode)
    assert width.cd_result.fully_cleared and dose.cd_result.fully_cleared
    assert width.cd_result.cd_um == pytest.approx(target, abs=1e-8)
    assert dose.cd_result.cd_um == pytest.approx(target, abs=1e-8)
    assert width.reinforcement_percent == pytest.approx(100*width.write_width_um/target)
    assert dose.reinforcement_percent == pytest.approx(100*dose.dose_scale)
