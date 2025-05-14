#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""This is to test the function _get_d_spacing in the interface.py file."""
import numpy as np
from packages.brillouin_calculator.core import _get_d_spacing


def test_cubic_spacing():
    a = b = c = 5.0
    alpha = beta = gamma = 90.0

    h, k, l = 1, 1, 1
    d_calculated = _get_d_spacing(h, k, l, a, b, c, alpha, beta, gamma)
    d_expected = a / np.sqrt(h**2 + k**2 + l**2)
    assert np.isclose(d_calculated, d_expected, atol=1e-5)

    h, k, l = 2, 0, 0
    d_calculated = _get_d_spacing(h, k, l, a, b, c, alpha, beta, gamma)
    d_expected = a / np.sqrt(h**2 + k**2 + l**2)
    assert np.isclose(d_calculated, d_expected, atol=1e-5)


def test_tetragonal_spacing():
    a = b = 4.0
    c = 6.0
    alpha = beta = gamma = 90.0

    h, k, l = 1, 1, 1
    d_calculated = _get_d_spacing(h, k, l, a, b, c, alpha, beta, gamma)
    d_expected = 1 / np.sqrt((h**2 + k**2) / a**2 + l**2 / c**2)
    assert np.isclose(d_calculated, d_expected, atol=1e-5)


def test_hexagonal_spacing():
    a = b = 3.0
    c = 5.0
    alpha = beta = 90.0
    gamma = 120.0

    h, k, l = 1, 1, 0
    d_calculated = _get_d_spacing(h, k, l, a, b, c, alpha, beta, gamma)
    d_expected = 1 / np.sqrt((4.0 / 3.0) * ((h**2 + h * k + k**2) / a**2) + l**2 / c**2)
    assert np.isclose(d_calculated, d_expected, atol=1e-5)


def test_orthorhombic_spacing():
    a = 3.0
    b = 4.0
    c = 5.0
    alpha = beta = gamma = 90.0

    h, k, l = 1, 1, 1
    d_calculated = _get_d_spacing(h, k, l, a, b, c, alpha, beta, gamma)
    d_expected = 1 / np.sqrt(h**2 / a**2 + k**2 / b**2 + l**2 / c**2)
    assert np.isclose(d_calculated, d_expected, atol=1e-5)


def test_monoclinic_spacing():
    a = 3.0
    b = 4.0
    c = 5.0
    alpha = 90.0
    beta = 100.0
    gamma = 90.0

    h, k, l = 1, 1, 1
    d_calculated = _get_d_spacing(h, k, l, a, b, c, alpha, beta, gamma)
    beta_rad = np.radians(beta)
    d_expected = 1 / np.sqrt(
        (1 / np.sin(beta_rad) ** 2)
        * (
            h**2 / a**2
            + (k**2 * np.sin(beta_rad) ** 2) / b**2
            + l**2 / c**2
            - (2 * h * l * np.cos(beta_rad)) / (a * c)
        )
    )
    assert np.isclose(d_calculated, d_expected, atol=1e-5)


def test_rhombohedral_spacing():
    a = b = c = 5.0
    alpha = beta = gamma = 70.0

    h, k, l = 1, 1, 1
    d_calculated = _get_d_spacing(h, k, l, a, b, c, alpha, beta, gamma)
    alpha_rad = np.radians(alpha)
    term1 = (h**2 + k**2 + l**2) * np.sin(alpha_rad) ** 2
    term2 = 2 * (h * k + k * l + h * l) * (np.cos(alpha_rad) ** 2 - np.cos(alpha_rad))
    denominator = a**2 * (1 - 3 * np.cos(alpha_rad) ** 2 + 2 * np.cos(alpha_rad) ** 3)
    d_expected = 1 / np.sqrt((term1 + term2) / denominator)
    assert np.isclose(d_calculated, d_expected, atol=1e-5)


def test_triclinic_spacing():
    a = 3.0
    b = 4.0
    c = 5.0
    alpha = 80.0
    beta = 85.0
    gamma = 100.0

    h, k, l = 1, 1, 1
    d_calculated = _get_d_spacing(h, k, l, a, b, c, alpha, beta, gamma)

    alpha_rad = np.radians(alpha)
    beta_rad = np.radians(beta)
    gamma_rad = np.radians(gamma)

    cos_alpha = np.cos(alpha_rad)
    cos_beta = np.cos(beta_rad)
    cos_gamma = np.cos(gamma_rad)
    V = (
        a
        * b
        * c
        * np.sqrt(
            1
            - cos_alpha**2
            - cos_beta**2
            - cos_gamma**2
            + 2 * cos_alpha * cos_beta * cos_gamma
        )
    )

    S11 = b**2 * c**2 * np.sin(alpha_rad) ** 2
    S22 = a**2 * c**2 * np.sin(beta_rad) ** 2
    S33 = a**2 * b**2 * np.sin(gamma_rad) ** 2
    S12 = a * b * c**2 * (cos_alpha * cos_beta - cos_gamma)
    S23 = a**2 * b * c * (cos_beta * cos_gamma - cos_alpha)
    S13 = a * b**2 * c * (cos_gamma * cos_alpha - cos_beta)

    one_over_d2 = (1 / V**2) * (
        S11 * h**2
        + S22 * k**2
        + S33 * l**2
        + 2 * S12 * h * k
        + 2 * S23 * k * l
        + 2 * S13 * h * l
    )
    d_expected = 1 / np.sqrt(one_over_d2)

    assert np.isclose(d_calculated, d_expected, atol=1e-5)
