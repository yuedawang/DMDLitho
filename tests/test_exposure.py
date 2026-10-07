import numpy as np
import pytest
from scipy.integrate import quad

from dmdlitho import normalized_dose_profile, normalized_peak_dose


def test_symmetry_and_scalar_vector_agreement():
    positions = np.linspace(-15, 15, 91)
    dose = normalized_dose_profile(positions, 6, 2)
    np.testing.assert_allclose(dose, normalized_dose_profile(-positions, 6, 2), atol=1e-15)
    assert normalized_dose_profile(positions[9], 6, 2) == pytest.approx(dose[9])


@pytest.mark.parametrize("width", [0.1, 1, 5.4, 6.6, 15, 40])
def test_equations_4_and_7(width):
    assert normalized_dose_profile(0, width, 2) == pytest.approx(normalized_peak_dose(width, 2))


def test_profile_matches_independent_gaussian_integral():
    for x in [0, 2, 4, 8]:
        expected = quad(lambda y: np.exp(-((x-y)/2)**2/2) / (2*np.sqrt(2*np.pi)), -3, 3)[0]
        assert normalized_dose_profile(x, 6, 2) == pytest.approx(expected, abs=1e-13)


def test_wide_feature_and_tail():
    assert normalized_peak_dose(40, 2) == pytest.approx(1, abs=1e-12)
    assert normalized_dose_profile(25, 6, 2) > 0  # Avoid subtracting two rounded erfs.


@pytest.mark.parametrize("bad", [0, -1, np.nan, np.inf])
def test_invalid_lengths(bad):
    with pytest.raises(ValueError):
        normalized_dose_profile(0, bad, 2)
    with pytest.raises(ValueError):
        normalized_peak_dose(6, bad)


def test_invalid_positions_and_peak_arrays():
    with pytest.raises(ValueError):
        normalized_dose_profile([0, np.nan], 6, 2)
    with pytest.raises(ValueError):
        normalized_peak_dose([1, 0], 2)
