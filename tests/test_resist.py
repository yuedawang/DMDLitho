import numpy as np
import pytest

from dmdlitho import (LithographyParameters, literal_thresholds,
                      reconstructed_thresholds, normalized_removed_depth)


@pytest.mark.parametrize("gamma", [0.46, 0.47, 0.48, 0.49, 0.5])
def test_threshold_normalizations(gamma):
    threshold, clearing = reconstructed_thresholds(gamma)
    assert threshold + clearing == pytest.approx(1)
    assert 1 / np.log(clearing / threshold) == pytest.approx(gamma)
    literal_t, literal_c = literal_thresholds(gamma)
    assert literal_t == 0.1
    assert literal_c == pytest.approx(0.1 * np.exp(1/gamma))
    assert literal_t != threshold


def test_removed_depth_convention_and_intermediate_dose():
    threshold, clearing = reconstructed_thresholds(0.46)
    intermediate = np.sqrt(threshold * clearing)
    doses = [0, threshold/2, threshold, intermediate, clearing, 2*clearing]
    np.testing.assert_allclose(normalized_removed_depth(doses, threshold, clearing, 0.46),
                               [0, 0, 0, 0.5, 1, 1], atol=1e-14)
    assert normalized_removed_depth(intermediate, threshold, clearing, 0.46) == pytest.approx(0.5)


@pytest.mark.parametrize("bad", [-1, np.nan, np.inf])
def test_invalid_dose(bad):
    threshold, clearing = reconstructed_thresholds(0.46)
    with pytest.raises(ValueError):
        normalized_removed_depth(bad, threshold, clearing, 0.46)


@pytest.mark.parametrize("bad", [0, -1, np.nan, np.inf, 1e-6, 1e20])
def test_invalid_or_unrepresentable_gamma(bad):
    for function in [literal_thresholds, reconstructed_thresholds]:
        with pytest.raises(ValueError):
            function(bad)


def test_inconsistent_resist_inputs_and_parameters():
    with pytest.raises(ValueError, match="Eq. .5."):
        normalized_removed_depth(0.5, 0.1, 0.9, 0.5)
    with pytest.raises(ValueError):
        normalized_removed_depth(0.5, 0.9, 0.1, 0.5)
    for name in ["gamma", "gaussian_radius_um", "resist_thickness_um"]:
        with pytest.raises(ValueError):
            LithographyParameters(**{name: 0})
