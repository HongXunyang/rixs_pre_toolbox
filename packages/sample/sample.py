"""This is a module for the sample class."""

import numpy as np


class Sample:
    """This is a class for the sample."""

    def __init__(self):
        """Initialize the sample."""
        self.parameters = None
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

    def initialize(self, parameters: dict, sample_to_lab_conversion):
        """Initialize the sample.

        sample_to_lab_conversion: function that converts vectors from sample coordinate system to
        lab coordinate system. This need to be included to the calculator class.
        """
        self.parameters = parameters

        # calculate vectors in sample coordinate system
        a, b, c, alpha, beta, gamma = self.get_lattice_parameters()
        self.a_vec_sample, self.b_vec_sample, self.c_vec_sample = (
            _get_real_space_vectors(a, b, c, alpha, beta, gamma)
        )
        self.a_star_vec_sample, self.b_star_vec_sample, self.c_star_vec_sample = (
            _get_reciprocal_space_vectors(a, b, c, alpha, beta, gamma)
        )
        self.a_vec_lab, self.b_vec_lab, self.c_vec_lab = sample_to_lab_conversion(
            self.a_vec_sample, self.b_vec_sample, self.c_vec_sample
        )
        self.a_star_vec_lab, self.b_star_vec_lab, self.c_star_vec_lab = (
            sample_to_lab_conversion(
                self.a_star_vec_sample, self.b_star_vec_sample, self.c_star_vec_sample
            )
        )

    def get_lattice_parameters(self):
        """Get the parameters of the sample. if None, raise "initialize the sample first"."""
        try:
            a, b, c = self.parameters["a"], self.parameters["b"], self.parameters["c"]
            alpha, beta, gamma = (
                self.parameters["alpha"],
                self.parameters["beta"],
                self.parameters["gamma"],
            )
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
