"""This is to test the function calculate_hkl in the interface.py file."""

import numpy as np

import pytest
from packages.brillouin_calculator.interface import BrillouinCalculator


@pytest.fixture
def calculator():
    """Returns a BrillouinCalculator instance (not yet initialized)"""
    return BrillouinCalculator()


# Orthorhombic case
params_orthorhombic = {
    "energy": 950,
    "a": 4.0,
    "b": 5.0,
    "c": 12.0,
    "alpha": 90,
    "beta": 90,
    "gamma": 90,
    "yaw": 0,
    "pitch": 0,
    "roll": 0,
    "theta": 0,
    "phi": 0,
    "chi": 0,
}


@pytest.mark.parametrize(
    "tth, theta, phi, chi, H_expected, K_expected, L_expected",
    [
        (150, 0, 0, 0, -0.5719, 0.0, -0.4597),
        (150, 90, 0, 0, 0.1532, 0.0, -1.716),
        (150, 90, 30, 23, 0.02098, -0.3377, -1.579),
    ],
)
def test_orthorhombic_calculate_hkl(
    calculator, tth, theta, phi, chi, H_expected, K_expected, L_expected
):
    """test the function calculate_hkl from the brillouin_calculator/interface.py in the orthorhombic case"""
    calculator.initialize(params_orthorhombic)
    results = calculator.calculate_hkl(tth, theta, phi, chi)
    assert np.isclose(results["H"], H_expected, atol=1e-3)
    assert np.isclose(results["K"], K_expected, atol=1e-3)
    assert np.isclose(results["L"], L_expected, atol=1e-3)


params_monoclinic = {
    "energy": 950,
    "a": 4.0,
    "b": 5.0,
    "c": 12.0,
    "alpha": 90,
    "beta": 90,
    "gamma": 120,
    "yaw": 0,
    "pitch": 0,
    "roll": 0,
    "theta": 0,
    "phi": 0,
    "chi": 0,
}


@pytest.mark.parametrize(
    "tth, theta, phi, chi, H_expected, K_expected, L_expected",
    [
        (150, 0, 0, 0, -0.5719, 0.0, -0.4597),
    ],
)
def test_monoclinic_calculate_hkl(
    calculator, tth, theta, phi, chi, H_expected, K_expected, L_expected
):
    """test the function calculate_hkl from the brillouin_calculator/interface.py in the monoclinic case"""
    pass
