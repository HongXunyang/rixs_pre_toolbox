"""This is to test the function calculate_angles in the interface.py file."""

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
    "H, K, L, fixed_angle, fixed_angle_name, tth_expected, theta_expected, phi_expected, chi_expected",
    [
        (-0.5719, 0.0, -0.4597, 0, "chi", 150.0, 0.0, 0.0, 0.0),
        (0.1532, 0.0, -1.716, 0, "chi", 150.0, 90.0, 0.0, 0.0),
        (0.02098, -0.3377, -1.579, 23, "chi", 150.0, 90.0, 30.0, 23),
    ],
)
def test_orthorhombic_calculate_angles(
    calculator,
    H,
    K,
    L,
    fixed_angle,
    fixed_angle_name,
    tth_expected,
    theta_expected,
    phi_expected,
    chi_expected,
):
    """test the function calculate_angles from the brillouin_calculator/interface.py in the orthorhombic case"""
    calculator.initialize(params_orthorhombic)
    results = calculator.calculate_angles(H, K, L, fixed_angle, fixed_angle_name)
    print(
        f"tth: {results['tth']}, tth_expected: {tth_expected},\n theta: {results['theta']}, theta_expected: {theta_expected},\n phi: {results['phi']}, phi_expected: {phi_expected},\n chi: {results['chi']}, chi_expected: {chi_expected}"
    )
    assert np.isclose(results["tth"], tth_expected, atol=1.0)
    assert np.isclose(results["theta"], theta_expected, atol=1.0)
    assert np.isclose(results["phi"], phi_expected, atol=1.0)
    assert np.isclose(results["chi"], chi_expected, atol=1.0)
