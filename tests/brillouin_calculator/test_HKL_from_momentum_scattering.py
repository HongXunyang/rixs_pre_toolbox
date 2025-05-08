"""This is to test the function _get_HKL_from_momentum_scattering in the interface.py file."""

import numpy as np
import pytest
from packages.brillouin_calculator.interface import (
    _get_HKL_from_momentum_scattering,
    BrillouinCalculator,
)


@pytest.fixture
def calculator():
    # Returns a BrillouinCalculator instance (not yet initialized)
    return BrillouinCalculator()


e_H = np.array([1, 0, 0])
e_K = np.array([0, 1, 0])
e_L = np.array([0, 0, 1])


@pytest.mark.parametrize(
    "params, tth, theta, phi, chi, H_expected, K_expected, L_expected",
    [
        (
            {
                "a": 5,
                "b": 4,
                "c": 12,
                "alpha": 90,
                "beta": 90,
                "gamma": 90,
                "energy": 950,
                "e_H": e_H,
                "e_K": e_K,
                "e_L": e_L,
            },
            152,
            80,
            0,
            0,
            0.0518615991,
            0,
            -1.779973,
        ),
        # Add more parameter sets here
    ],
)
def test_HKL_from_momentum_scattering(
    calculator, params, tth, theta, phi, chi, H_expected, K_expected, L_expected
):
    calculator.initialize(params)

    # test the function _get_HKL_from_momentum_scattering
    momentum = np.array(
        [
            2 * np.pi * H_expected / (params["a"]),
            2 * np.pi * K_expected / (params["b"]),
            2 * np.pi * L_expected / (params["c"]),
        ]
    )
    print(f"momentum: {momentum}")
    a_vec, b_vec, c_vec = calculator.a_vec, calculator.b_vec, calculator.c_vec
    H_test, K_test, L_test = _get_HKL_from_momentum_scattering(
        momentum, a_vec, b_vec, c_vec
    )

    assert np.isclose(H_expected, H_test)
    assert np.isclose(K_expected, K_test)
    assert np.isclose(L_expected, L_test)
