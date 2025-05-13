"""
This is a test file to test the methods `get_lattice_basis` and `get_sample_basis` in the `Lab` and
the `Sample` class.
"""

import pytest
import numpy as np
import sys
import os

sys.path.append(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
)
from packages.classes.lab import Lab
from packages.classes.sample import Sample
from packages.utils.translation import angle_to_matrix, euler_to_matrix


params = {
    "a": 1,
    "b": 3,
    "c": 5,
    "alpha": 95,
    "beta": 120,
    "gamma": 102,
    "roll": 0,
    "pitch": 90,
    "yaw": 0,
    "theta": 90,
    "phi": 0,
    "chi": 0,
}

params_2 = {
    "a": 1,
    "b": 3,
    "c": 5,
    "alpha": 95,
    "beta": 120,
    "gamma": 102,
    "roll": 90,
    "pitch": 0,
    "yaw": 90,
    "theta": 90,
    "phi": 90,
    "chi": 0,
}


@pytest.fixture
def lab():
    return Lab()


def test_sample_and_lattice_basis_orthogonality(lab):
    """Test that the lab basis vectors are orthogonal."""
    lab.initialize(**params)
    ex_sample, ey_sample, ez_sample = lab.get_sample_basis()
    ex_lattice, ey_lattice, ez_lattice = lab.get_lattice_basis()

    # Test orthogonality
    assert np.isclose(np.dot(ex_sample, ey_sample), 0)
    assert np.isclose(np.dot(ex_sample, ez_sample), 0)
    assert np.isclose(np.dot(ey_sample, ez_sample), 0)

    # Test normalization
    assert np.isclose(np.linalg.norm(ex_sample), 1)
    assert np.isclose(np.linalg.norm(ey_sample), 1)
    assert np.isclose(np.linalg.norm(ez_sample), 1)

    # Test right-handedness
    assert np.allclose(np.cross(ex_sample, ey_sample), ez_sample)

    # Test orthogonality
    assert np.isclose(np.dot(ex_lattice, ey_lattice), 0)
    assert np.isclose(np.dot(ex_lattice, ez_lattice), 0)
    assert np.isclose(np.dot(ey_lattice, ez_lattice), 0)

    # Test normalization
    assert np.isclose(np.linalg.norm(ex_lattice), 1)
    assert np.isclose(np.linalg.norm(ey_lattice), 1)
    assert np.isclose(np.linalg.norm(ez_lattice), 1)

    # Test right-handedness
    assert np.allclose(np.cross(ex_lattice, ey_lattice), ez_lattice)


def test_lab_sample_basis_rotation(lab):
    """Test that the sample basis in lab frame is correctly rotated."""

    # test 1:
    lab.initialize(**params)

    ex_sample_in_lab, ey_sample_in_lab, ez_sample_in_lab = lab.get_sample_basis()

    # The sample basis in the lab frame should be the lab basis rotated by the rotation matrix
    expected_ex = np.array([0, 0, -1])
    expected_ey = np.array([0, 1, 0])
    expected_ez = np.array([1, 0, 0])

    assert np.allclose(ex_sample_in_lab, expected_ex)
    assert np.allclose(ey_sample_in_lab, expected_ey)
    assert np.allclose(ez_sample_in_lab, expected_ez)

    # test 2: roll = yaw = theta = phi = 90
    lab.initialize(**params_2)
    ex_sample_in_lab, ey_sample_in_lab, ez_sample_in_lab = lab.get_sample_basis()

    expected_ex = np.array([0, 1, 0])
    expected_ey = np.array([0, 0, 1])
    expected_ez = np.array([1, 0, 0])

    assert np.allclose(ex_sample_in_lab, expected_ex)
    assert np.allclose(ey_sample_in_lab, expected_ey)
    assert np.allclose(ez_sample_in_lab, expected_ez)


def test_sample_lattice_basis_rotation(lab):
    """Test that the lattice basis in sample frame is correctly rotated."""

    # test 1: theta = pitch = 90
    lab.initialize(**params)
    sample = lab.sample
    ex_lattice, ey_lattice, ez_lattice = sample.get_lattice_basis()

    # The lattice basis in the sample frame should be the sample basis rotated by the rotation matrix
    expected_ex = np.array([0, 0, -1])
    expected_ey = np.array([0, 1, 0])
    expected_ez = np.array([1, 0, 0])

    assert np.allclose(ex_lattice, expected_ex)
    assert np.allclose(ey_lattice, expected_ey)
    assert np.allclose(ez_lattice, expected_ez)

    # test 2: roll = yaw = theta = phi = 90
    lab.initialize(**params_2)
    sample = lab.sample
    ex_lattice, ey_lattice, ez_lattice = sample.get_lattice_basis()

    expected_ex = np.array([0, 1, 0])
    expected_ey = np.array([0, 0, 1])
    expected_ez = np.array([1, 0, 0])

    assert np.allclose(ex_lattice, expected_ex)
    assert np.allclose(ey_lattice, expected_ey)
    assert np.allclose(ez_lattice, expected_ez)


def test_lab_lattice_basis_rotation(lab):
    """Test that the lattice basis in lab frame is correctly rotated."""
    lab.initialize(**params)

    ex_lattice_in_lab, ey_lattice_in_lab, ez_lattice_in_lab = lab.get_lattice_basis()

    # The lattice basis in the lab frame should be the lattice basis in the sample frame rotated by the lab rotation matrix
    expected_ex = np.array([-1, 0, 0])
    expected_ey = np.array([0, 1, 0])
    expected_ez = np.array([0, 0, -1])

    assert np.allclose(ex_lattice_in_lab, expected_ex)
    assert np.allclose(ey_lattice_in_lab, expected_ey)
    assert np.allclose(ez_lattice_in_lab, expected_ez)

    # test 2: roll = yaw = theta = phi = 90
    lab.initialize(**params_2)
    ex_lattice, ey_lattice, ez_lattice = lab.get_lattice_basis()

    expected_ex = np.array([0, 0, 1])
    expected_ey = np.array([1, 0, 0])
    expected_ez = np.array([0, 1, 0])

    assert np.allclose(ex_lattice, expected_ex)
    assert np.allclose(ey_lattice, expected_ey)
    assert np.allclose(ez_lattice, expected_ez)
