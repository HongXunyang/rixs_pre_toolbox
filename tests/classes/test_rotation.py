"""Tests for the sample and lab rotation functionality."""

import numpy as np
import sys
import os

# Add the parent directory to sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from packages.classes.sample import Sample
from packages.classes.lab import Lab
from packages.utils import euler_to_matrix, angle_to_matrix


def test_sample_initialization():
    """Test if Sample initializes correctly."""
    sample = Sample()
    sample.initialize(1.0, 2.0, 3.0, 90.0, 90.0, 90.0, 0.0, 0.0, 0.0)

    a, b, c, alpha, beta, gamma = sample.get_lattice_parameters()
    assert a == 1.0
    assert b == 2.0
    assert c == 3.0
    assert alpha == 90.0
    assert beta == 90.0
    assert gamma == 90.0

    assert sample.roll == 0.0
    assert sample.pitch == 0.0
    assert sample.yaw == 0.0


def test_lab_initialization():
    """Test if Lab initializes correctly."""
    lab = Lab()
    lab.initialize(1.0, 2.0, 3.0, 90.0, 90.0, 90.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    a, b, c, alpha, beta, gamma = lab.get_lattice_parameters()
    assert a == 1.0
    assert b == 2.0
    assert c == 3.0
    assert alpha == 90.0
    assert beta == 90.0
    assert gamma == 90.0

    assert lab.theta == 0.0
    assert lab.phi == 0.0
    assert lab.chi == 0.0


def test_sample_orthorhombic_no_rotation():
    """Test orthorhombic sample with no rotation."""
    sample = Sample()
    sample.initialize(1.0, 2.0, 3.0, 90.0, 90.0, 90.0, 0.0, 0.0, 0.0)

    a_vec, b_vec, c_vec = sample.get_real_space_vectors()

    # For orthorhombic cell with no rotation, vectors should align with axes
    np.testing.assert_almost_equal(a_vec, np.array([1.0, 0.0, 0.0]))
    np.testing.assert_almost_equal(b_vec, np.array([0.0, 2.0, 0.0]))
    np.testing.assert_almost_equal(c_vec, np.array([0.0, 0.0, 3.0]))

    # Check reciprocal vectors
    a_star_vec, b_star_vec, c_star_vec = sample.get_reciprocal_space_vectors()
    np.testing.assert_almost_equal(a_star_vec, np.array([1.0, 0.0, 0.0]) * np.pi * 2)
    np.testing.assert_almost_equal(b_star_vec, np.array([0.0, 0.5, 0.0]) * np.pi * 2)
    np.testing.assert_almost_equal(
        c_star_vec, np.array([0.0, 0.0, 1.0 / 3.0]) * np.pi * 2
    )


def test_sample_with_rotation():
    """Test sample with rotation."""
    sample = Sample()
    # Orthorhombic cell with 90 degree rotation around z-axis (yaw)
    sample.initialize(1.0, 2.0, 3.0, 90.0, 90.0, 90.0, 0.0, 0.0, 90.0)

    a_vec, b_vec, c_vec = sample.get_real_space_vectors()

    # After 90 degree rotation around z, a should point along -y and b along x
    np.testing.assert_almost_equal(a_vec, np.array([0.0, 1.0, 0.0]))
    np.testing.assert_almost_equal(b_vec, np.array([-2.0, 0.0, 0.0]))
    np.testing.assert_almost_equal(c_vec, np.array([0.0, 0.0, 3.0]))


def test_lab_with_rotation():
    """Test lab with rotation."""
    lab = Lab()
    # Orthorhombic cell, no sample rotation, but with lab rotation in y-axis (theta)
    lab.initialize(1.0, 2.0, 3.0, 90.0, 90.0, 90.0, 0.0, 0.0, 0.0, -90, 0.0, 0)

    a_vec_lab, b_vec_lab, c_vec_lab = lab.get_real_space_vectors()

    # After 90 degree chi rotation, a should point along -y and b along x
    np.testing.assert_almost_equal(a_vec_lab, np.array([0.0, 0.0, 1.0]))
    np.testing.assert_almost_equal(b_vec_lab, np.array([0.0, 2, 0.0]))
    np.testing.assert_almost_equal(c_vec_lab, np.array([-3.0, 0.0, 0.0]))


def test_combined_rotations():
    """Test combined sample and lab rotations."""
    lab = Lab()
    # Orthorhombic cell with sample and lab rotations, sample rotation in y-axis (pitch) and lab rotation in x-axis (chi)
    lab.initialize(1.0, 2.0, 3.0, 90.0, 90.0, 90.0, 0.0, 90.0, 0.0, 0.0, 0.0, 90.0)

    a_vec_lab, b_vec_lab, c_vec_lab = lab.get_real_space_vectors()

    # Sample pitch=90 rotates a to -z, c to x
    # Then lab theta=90 rotates these again
    np.testing.assert_almost_equal(a_vec_lab, np.array([0.0, 1.0, 0.0]))
    np.testing.assert_almost_equal(b_vec_lab, np.array([0.0, 0.0, 2.0]))
    np.testing.assert_almost_equal(c_vec_lab, np.array([3.0, 0.0, 0.0]))


def test_rotation_matrix_consistency():
    """Test that our rotation matrices are consistent with vector transformation."""
    # Test for euler_to_matrix
    roll, pitch, yaw = 30.0, 45.0, 60.0
    rotation_matrix = euler_to_matrix(roll, pitch, yaw)

    # Create a test vector
    vec = np.array([1.0, 0.0, 0.0])
    rotated_vec = rotation_matrix @ vec

    # Create a sample and apply the same rotation
    sample = Sample()
    sample.initialize(1.0, 1.0, 1.0, 90.0, 90.0, 90.0, roll, pitch, yaw)
    a_vec, _, _ = sample.get_real_space_vectors()

    np.testing.assert_almost_equal(rotated_vec, a_vec)

    # Test for angle_to_matrix
    theta, phi, chi = -30.0, 45.0, 60.0
    lab_rotation = angle_to_matrix(theta, phi, chi)

    lab_rotated_vec = lab_rotation @ vec

    # Create a lab and apply the same rotation
    lab = Lab()
    lab.initialize(1.0, 1.0, 1.0, 90.0, 90.0, 90.0, 0.0, 0.0, 0.0, theta, phi, chi)
    a_vec_lab, _, _ = lab.get_real_space_vectors()

    np.testing.assert_almost_equal(lab_rotated_vec, a_vec_lab)


def test_non_orthogonal_lattice():
    """Test rotation with a non-orthogonal lattice."""
    # Create a monoclinic cell
    sample = Sample()
    sample.initialize(1.0, 2.0, 3.0, 90.0, 90.0, 120.0, 0.0, 0.0, 0.0)

    a_vec, b_vec, c_vec = sample.get_real_space_vectors()

    # Monoclinic cell should have c with an x component
    np.testing.assert_almost_equal(a_vec, np.array([1.0, 0.0, 0.0]))
    np.testing.assert_almost_equal(
        b_vec,
        np.array([2 * np.cos(np.radians(120)), 2.0 * np.sin(np.radians(120)), 0.0]),
    )
    np.testing.assert_almost_equal(c_vec, np.array([0.0, 0.0, 3.0]))

    # Now rotate this cell about the x-axis by 90 degrees
    sample.initialize(1.0, 2.0, 3.0, 90.0, 90.0, 120.0, 90.0, 0.0, 0)
    a_vec, b_vec, c_vec = sample.get_real_space_vectors()

    np.testing.assert_almost_equal(a_vec, np.array([1, 0, 0.0]))
    np.testing.assert_almost_equal(
        b_vec,
        np.array([2 * np.cos(np.radians(120)), 0, 2.0 * np.sin(np.radians(120))]),
    )
    np.testing.assert_almost_equal(c_vec, np.array([0, -3, 0]))


def test_reciprocal_space_rotation():
    """Test that reciprocal space vectors rotate correctly."""
    lab = Lab()
    lab.initialize(1.0, 2.0, 3.0, 90.0, 90.0, 90.0, 0.0, 0.0, 0.0, 0.0, 0.0, 90.0)

    a_star_vec, b_star_vec, c_star_vec = lab.get_reciprocal_space_vectors()

    np.testing.assert_almost_equal(a_star_vec, np.array([1.0, 0.0, 0.0]) * np.pi * 2)
    np.testing.assert_almost_equal(b_star_vec, np.array([0, 0, 0.5]) * np.pi * 2)
    np.testing.assert_almost_equal(c_star_vec, np.array([0.0, -1 / 3, 0]) * np.pi * 2)
