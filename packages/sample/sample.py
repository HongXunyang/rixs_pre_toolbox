"""This is a module for the sample class."""

import numpy as np


class Sample:
    """This is a class for the sample."""

    def __init__(self):
        """Initialize the sample."""
        self.a = None
        self.b = None
        self.c = None
        self.alpha = None
        self.beta = None
        self.gamma = None
        # vectors in sample coordinate system
        self.a_vec_sample = None
        self.b_vec_sample = None
        self.c_vec_sample = None
        self.a_star_vec_sample = None
        self.b_star_vec_sample = None
        self.c_star_vec_sample = None

        # vectors in lab coordinate system
        self.a_vec_lab = None
        self.b_vec_lab = None
        self.c_vec_lab = None
        self.a_star_vec_lab = None
        self.b_star_vec_lab = None
        self.c_star_vec_lab = None

    def initialize(self, a, b, c, alpha, beta, gamma, roll, pitch, yaw):
        """Initialize the sample.

        Args:
            a, b, c (float): Lattice constants in Angstroms
            alpha, beta, gamma (float): Lattice angles in degrees
            roll, pitch, yaw (float): Euler angles in degrees
        """
        # First set the lattice parameters
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma

        # Then calculate vectors in sample coordinate system
        self.a_vec_sample, self.b_vec_sample, self.c_vec_sample = (
            _get_real_space_vectors(a, b, c, alpha, beta, gamma)
        )
        self.a_star_vec_sample, self.b_star_vec_sample, self.c_star_vec_sample = (
            _get_reciprocal_space_vectors(a, b, c, alpha, beta, gamma)
        )

        self.a_vec_lab, self.b_vec_lab, self.c_vec_lab = _sample_to_lab_conversion(
            self.a_vec_sample, self.b_vec_sample, self.c_vec_sample, roll, pitch, yaw
        )
        self.a_star_vec_lab, self.b_star_vec_lab, self.c_star_vec_lab = (
            _sample_to_lab_conversion(
                self.a_star_vec_sample,
                self.b_star_vec_sample,
                self.c_star_vec_sample,
                roll,
                pitch,
                yaw,
            )
        )

    def get_lattice_parameters(self):
        """Get the parameters of the sample. if None, raise "initialize the sample first"."""
        try:
            a, b, c = self.a, self.b, self.c
            alpha, beta, gamma = self.alpha, self.beta, self.gamma
            return a, b, c, alpha, beta, gamma
        except KeyError as exc:
            raise ValueError("initialize the sample first") from exc

    def get_real_space_vectors(self, frame: str = "sample"):
        """Get the real space vectors.
        Args:
            frame (str): "sample" or "lab"
        """
        if frame == "sample":
            return self.a_vec_sample, self.b_vec_sample, self.c_vec_sample
        elif frame == "lab":
            return self.a_vec_lab, self.b_vec_lab, self.c_vec_lab
        else:
            raise ValueError("invalid frame, try 'sample' or 'lab'")

    def get_reciprocal_space_vectors(self, frame: str = "sample"):
        """Get the reciprocal space vectors.
        Args:
            frame (str): "sample" or "lab"
        """
        if frame == "sample":
            return (
                self.a_star_vec_sample,
                self.b_star_vec_sample,
                self.c_star_vec_sample,
            )
        elif frame == "lab":
            return self.a_star_vec_lab, self.b_star_vec_lab, self.c_star_vec_lab
        else:
            raise ValueError("invalid frame, try 'sample' or 'lab'")


def _get_reciprocal_space_vectors(a, b, c, alpha, beta, gamma):
    """Get the reciprocal space vectors a_star_vec, b_star_vec, c_star_vec from the lattice
    parameters, angles in degrees. These vectors are in the crystal coordinate system.
    """
    a_vec, b_vec, c_vec = _get_real_space_vectors(a, b, c, alpha, beta, gamma)
    volumn = abs(np.dot(a_vec, np.cross(b_vec, c_vec)))
    a_star_vec = 2 * np.pi * np.cross(b_vec, c_vec) / volumn
    b_star_vec = 2 * np.pi * np.cross(c_vec, a_vec) / volumn
    c_star_vec = 2 * np.pi * np.cross(a_vec, b_vec) / volumn
    return a_star_vec, b_star_vec, c_star_vec


def _get_real_space_vectors(a, b, c, alpha, beta, gamma):
    """Get the real space vectors a_vec, b_vec, c_vec from the lattice parameters.
    - a_vec is by-default along x-axis (a, 0, 0)
    - b_vec is by-default (b cos gamma, b sin gamma, 0) on the x-y plane,
    - c_vec is then calculated
    The above convention defines the crystal coordinate system.

    Args:
        a, b, c (float): Lattice constants in Angstroms
        alpha, beta, gamma (float): Lattice angles in degrees

    Returns:
        a_vec, b_vec, c_vec (np.ndarray): Real space vectors
    """
    alpha_rad, beta_rad, gamma_rad = (
        np.radians(alpha),
        np.radians(beta),
        np.radians(gamma),
    )
    a_vec = np.array([a, 0, 0])
    b_vec = np.array([b * np.cos(gamma_rad), b * np.sin(gamma_rad), 0])
    c_vec_x = c * np.cos(beta_rad)
    c_vec_y = (
        c
        * (np.cos(alpha_rad) - np.cos(beta_rad) * np.cos(gamma_rad))
        / np.sin(gamma_rad)
    )
    c_vec_z = np.sqrt(c**2 - c_vec_x**2 - c_vec_y**2)
    c_vec = np.array([c_vec_x, c_vec_y, c_vec_z])
    return a_vec, b_vec, c_vec


def euler_to_matrix(roll, pitch, yaw):
    """Convert Euler angles to rotation matrix. We follows the ZYX convention.

    Args:
        roll (float): rotation about the new X axis in degrees
        pitch (float): rotation about the new Y axis in degrees
        yaw (float): rotation about the original z axis in degrees

    Returns:
        rotation_matrix (np.ndarray): Rotation matrix
    """
    roll_rad, pitch_rad, yaw_rad = (
        np.radians(roll),
        np.radians(pitch),
        np.radians(yaw),
    )
    Rx = np.array(
        [
            [1, 0, 0],
            [0, np.cos(roll_rad), -np.sin(roll_rad)],
            [0, np.sin(roll_rad), np.cos(roll_rad)],
        ]
    )

    Ry = np.array(
        [
            [np.cos(pitch_rad), 0, np.sin(pitch_rad)],
            [0, 1, 0],
            [-np.sin(pitch_rad), 0, np.cos(pitch_rad)],
        ]
    )

    Rz = np.array(
        [
            [np.cos(yaw_rad), -np.sin(yaw_rad), 0],
            [np.sin(yaw_rad), np.cos(yaw_rad), 0],
            [0, 0, 1],
        ]
    )

    return Rz @ Ry @ Rx  # ZYX order


def _sample_to_lab_conversion(
    a_vec_sample, b_vec_sample, c_vec_sample, roll, pitch, yaw
):
    """Convert vectors from sample coordinate system to lab coordinate system.

    Args:
        a_vec_sample, b_vec_sample, c_vec_sample (np.ndarray): Vectors in sample coordinate system
        roll, pitch, yaw (float): Euler angles in degrees

    Returns:
        a_vec_lab, b_vec_lab, c_vec_lab (np.ndarray): Vectors in lab coordinate system
    """
    rotation_matrix = euler_to_matrix(roll, pitch, yaw)
    a_vec_lab = rotation_matrix @ a_vec_sample
    b_vec_lab = rotation_matrix @ b_vec_sample
    c_vec_lab = rotation_matrix @ c_vec_sample
    return a_vec_lab, b_vec_lab, c_vec_lab


def _lab_to_sample_conversion(a_vec_lab, b_vec_lab, c_vec_lab, roll, pitch, yaw):
    """Convert vectors from lab coordinate system to sample coordinate system.

    Args:
        a_vec_lab, b_vec_lab, c_vec_lab (np.ndarray): Vectors in lab coordinate system
        roll, pitch, yaw (float): Euler angles in degrees

    Returns:
        a_vec_sample, b_vec_sample, c_vec_sample (np.ndarray): Vectors in sample coordinate system
    """
    rotation_matrix = euler_to_matrix(roll, pitch, yaw)
    a_vec_sample = rotation_matrix.T @ a_vec_lab
    b_vec_sample = rotation_matrix.T @ b_vec_lab
    c_vec_sample = rotation_matrix.T @ c_vec_lab
    return a_vec_sample, b_vec_sample, c_vec_sample
