"""This is to test the function calculate_angles with tth fixed in the interface.py file."""

import numpy as np
import pytest
from packages.brillouin_calculator.interface import BrillouinCalculator


# TEST CASE 1: WITH ZERO CHI AND PHI
@pytest.fixture
def calculator():
    """Create and initialize a BrillouinCalculator instance."""
    calc = BrillouinCalculator()
    params = {
        "a": 4,
        "b": 4,
        "c": 12,
        "alpha": 90.0,
        "beta": 90.0,
        "gamma": 90.0,
        "energy": 930,
        "roll": 0.0,
        "pitch": 0.0,
        "yaw": 0.0,
    }
    calc.initialize(params)
    return calc


def test_calculate_angles_tth_fixed(calculator):
    """Test calculate_angles with tth fixed."""
    # Test parameters
    tth = 132.0
    H = -0.2230
    K = 0
    L = None  # L will be calculated based on tth
    fixed_angle = 0.0
    fixed_angle_name = "chi"

    # Expected angles that should be in the results
    expected_tth = tth  # Same as the fixed angle
    expected_theta = 42  # To be filled by user
    expected_phi = 0  # To be filled by user
    expected_chi = 0  # To be filled by user
    expected_L = -1.5024  # To be filled by user

    # Call the calculate_angles_tth_fixed method
    result = calculator.calculate_angles_tth_fixed(
        tth=tth,
        H=H,
        K=K,
        L=L,
        fixed_angle_name=fixed_angle_name,
        fixed_angle=fixed_angle,
    )
    tth_result = result["tth"]
    theta_result = result["theta"]
    phi_result = result["phi"]
    chi_result = result["chi"]
    L_result = result["L"]

    # Check if expected values are in the results
    # Note: Skip checks for None values
    if expected_tth is not None:
        assert any(
            np.isclose(tth_val, expected_tth, atol=0.1) for tth_val in tth_result
        ), f"Expected tth {expected_tth}, got {tth_result}"

    if expected_theta is not None:
        assert any(
            np.isclose(theta, expected_theta, atol=0.1) for theta in theta_result
        ), f"Expected theta {expected_theta}, got {theta_result}"

    if expected_phi is not None:
        assert any(
            np.isclose(phi, expected_phi, atol=0.1) for phi in phi_result
        ), f"Expected phi {expected_phi}, got {phi_result}"

    if expected_chi is not None:
        assert any(
            np.isclose(chi, expected_chi, atol=0.1) for chi in chi_result
        ), f"Expected chi {expected_chi}, got {chi_result}"

    if expected_L is not None:
        assert np.isclose(
            L_result, expected_L, atol=0.01
        ), f"Expected L {expected_L}, got {L_result}"


# TEST CASE 2: WITH NON-ZERO CHI AND PHI
@pytest.fixture
def calculator_rotated():
    """Create and initialize a BrillouinCalculator instance with non-zero roll/pitch/yaw."""
    calc = BrillouinCalculator()
    params = {
        "a": 4,
        "b": 4,
        "c": 12,
        "alpha": 90.0,
        "beta": 90.0,
        "gamma": 90.0,
        "energy": 930,
        "roll": 0,
        "pitch": 0,
        "yaw": 0,
    }
    calc.initialize(params)
    return calc


def test_calculate_angles_tth_fixed_rotated(calculator_rotated):
    """Test calculate_angles with tth fixed for rotated crystal."""
    # Test parameters
    tth = 111.0
    H = -0.3849
    K = -0.0916
    L = None  # L will be calculated based on tth
    fixed_angle = 44
    fixed_angle_name = "chi"

    # Expected angles that should be in the results
    expected_tth = tth  # Same as the fixed angle
    expected_theta = 22  # To be filled by user
    expected_phi = 33  # To be filled by user
    expected_chi = fixed_angle  # To be filled by user
    expected_L = -0.8899  # To be filled by user

    # Call the calculate_angles_tth_fixed method
    result = calculator_rotated.calculate_angles_tth_fixed(
        tth=tth,
        H=H,
        K=K,
        L=L,
        fixed_angle_name=fixed_angle_name,
        fixed_angle=fixed_angle,
    )
    tth_result = result["tth"]
    theta_result = result["theta"]
    phi_result = result["phi"]
    chi_result = result["chi"]
    L_result = result["L"]

    # Check if expected values are in the results
    # Note: Skip checks for None values
    if expected_tth is not None:
        assert any(
            np.isclose(tth_val, expected_tth, atol=0.1) for tth_val in tth_result
        ), f"Expected tth {expected_tth}, got {tth_result}"

    if expected_theta is not None:
        assert any(
            np.isclose(theta, expected_theta, atol=0.1) for theta in theta_result
        ), f"Expected theta {expected_theta}, got {theta_result}"

    if expected_phi is not None:
        assert any(
            np.isclose(phi, expected_phi, atol=0.1) for phi in phi_result
        ), f"Expected phi {expected_phi}, got {phi_result}"

    if expected_chi is not None:
        assert any(
            np.isclose(chi, expected_chi, atol=0.1) for chi in chi_result
        ), f"Expected chi {expected_chi}, got {chi_result}"

    if expected_L is not None:
        assert np.isclose(
            L_result, expected_L, atol=0.01
        ), f"Expected L {expected_L}, got {L_result}"
