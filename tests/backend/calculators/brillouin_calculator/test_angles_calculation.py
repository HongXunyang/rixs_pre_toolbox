"""This is to test the function calculate_angles with chi and phi fixed in the interface.py file."""

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


def test_calculate_angles_chi_fixed(calculator):
    """Test calculate_angles with chi fixed."""
    # Test parameters
    H = -0.2230
    K = 0
    L = -1.5024
    fixed_angle = 0.0
    fixed_angle_name = "chi"

    # Expected angles that should be in the results
    # Replace these with your actual expected values
    expected_tth = 132  # Replace with your expected value
    expected_theta = 42  # Replace with your expected value
    expected_phi = 0  # Replace with your expected value
    expected_chi = fixed_angle

    # Call the calculate_angles method
    result = calculator.calculate_angles(H, K, L, fixed_angle, fixed_angle_name)

    # Check if the calculation was successful
    assert result["success"] is True
    assert result["error"] is None

    # Get the results
    tth_results = result["tth"]
    theta_results = result["theta"]
    phi_results = result["phi"]
    chi_results = result["chi"]

    # Check if expected values are in the results
    # Note: Skip checks for None values
    if expected_tth is not None:
        assert any(
            np.isclose(tth, expected_tth, atol=0.1) for tth in tth_results
        ), f"Expected tth {expected_tth}, got {tth_results}"

    if expected_theta is not None:
        assert any(
            np.isclose(theta, expected_theta, atol=0.1) for theta in theta_results
        ), f"Expected theta {expected_theta}, got {theta_results}"

    if expected_phi is not None:
        assert any(
            np.isclose(phi, expected_phi, atol=0.1) for phi in phi_results
        ), f"Expected phi {expected_phi}, got {phi_results}"

    if expected_chi is not None:
        assert all(
            np.isclose(chi, expected_chi, atol=0.1) for chi in chi_results
        ), f"Expected chi {expected_chi}, got {chi_results}"


def test_calculate_angles_phi_fixed(calculator):
    """Test calculate_angles with phi fixed."""
    # Test parameters
    H = -0.2230
    K = 0
    L = -1.5024
    fixed_angle = 0
    fixed_angle_name = "phi"

    # Expected angles that should be in the results
    # Replace these with your actual expected values
    expected_tth = 132  # Replace with your expected value
    expected_theta = 42  # Replace with your expected value
    expected_phi = fixed_angle
    expected_chi = 0  # Replace with your expected value

    # Call the calculate_angles method
    result = calculator.calculate_angles(H, K, L, fixed_angle, fixed_angle_name)

    # Get the results
    tth_results = result["tth"]
    theta_results = result["theta"]
    phi_results = result["phi"]
    chi_results = result["chi"]

    # Check if expected values are in the results
    # Note: Skip checks for None values
    if expected_tth is not None:
        assert any(
            np.isclose(tth, expected_tth, atol=0.1) for tth in tth_results
        ), f"Expected tth {expected_tth}, got {tth_results}"

    if expected_theta is not None:
        assert any(
            np.isclose(theta, expected_theta, atol=0.1) for theta in theta_results
        ), f"Expected theta {expected_theta}, got {theta_results}"

    if expected_phi is not None:
        assert all(
            np.isclose(phi, expected_phi, atol=0.1) for phi in phi_results
        ), f"Expected phi {expected_phi}, got {phi_results}"

    if expected_chi is not None:
        assert any(
            np.isclose(chi, expected_chi, atol=0.1) for chi in chi_results
        ), f"Expected chi {expected_chi}, got {chi_results}"


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


def test_calculate_angles_chi_fixed_rotated(calculator_rotated):
    """Test calculate_angles with chi fixed for rotated crystal."""
    # Test parameters
    H = -0.3849
    K = -0.0916
    L = -0.8899
    fixed_angle = 44
    fixed_angle_name = "chi"

    # Expected angles that should be in the results
    expected_tth = 111
    expected_theta = 22  # Different due to crystal rotation
    expected_phi = 33  # Different due to crystal rotation
    expected_chi = fixed_angle

    # Call the calculate_angles method
    result = calculator_rotated.calculate_angles(H, K, L, fixed_angle, fixed_angle_name)

    # Check if the calculation was successful
    assert result["success"] is True
    assert result["error"] is None

    # Get the results
    tth_results = result["tth"]
    theta_results = result["theta"]
    phi_results = result["phi"]
    chi_results = result["chi"]

    # Check if expected values are in the results
    if expected_tth is not None:
        assert any(
            np.isclose(tth, expected_tth, atol=0.1) for tth in tth_results
        ), f"Expected tth {expected_tth}, got {tth_results}"

    if expected_theta is not None:
        assert any(
            np.isclose(theta, expected_theta, atol=0.1) for theta in theta_results
        ), f"Expected theta {expected_theta}, got {theta_results}"

    if expected_phi is not None:
        assert any(
            np.isclose(phi, expected_phi, atol=0.1) for phi in phi_results
        ), f"Expected phi {expected_phi}, got {phi_results}"

    if expected_chi is not None:
        assert all(
            np.isclose(chi, expected_chi, atol=0.1) for chi in chi_results
        ), f"Expected chi {expected_chi}, got {chi_results}"


def test_calculate_angles_phi_fixed_rotated(calculator_rotated):
    """Test calculate_angles with phi fixed for rotated crystal."""
    # Test parameters
    H = -0.3849
    K = -0.0916
    L = -0.8899
    fixed_angle = 33
    fixed_angle_name = "phi"

    # Expected angles that should be in the results
    expected_tth = 111
    expected_theta = 22  # Different due to crystal rotation
    expected_phi = fixed_angle
    expected_chi = 44  # Different due to crystal rotation

    # Call the calculate_angles method
    result = calculator_rotated.calculate_angles(H, K, L, fixed_angle, fixed_angle_name)

    # Get the results
    tth_results = result["tth"]
    theta_results = result["theta"]
    phi_results = result["phi"]
    chi_results = result["chi"]

    # Check if expected values are in the results
    if expected_tth is not None:
        assert any(
            np.isclose(tth, expected_tth, atol=0.1) for tth in tth_results
        ), f"Expected tth {expected_tth}, got {tth_results}"

    if expected_theta is not None:
        assert any(
            np.isclose(theta, expected_theta, atol=0.1) for theta in theta_results
        ), f"Expected theta {expected_theta}, got {theta_results}"

    if expected_phi is not None:
        assert all(
            np.isclose(phi, expected_phi, atol=0.1) for phi in phi_results
        ), f"Expected phi {expected_phi}, got {phi_results}"

    if expected_chi is not None:
        assert any(
            np.isclose(chi, expected_chi, atol=0.1) for chi in chi_results
        ), f"Expected chi {expected_chi}, got {chi_results}"
