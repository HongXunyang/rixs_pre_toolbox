"""This is a module for the sample class."""
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from packages.utils import (
    get_real_space_vectors,
    get_reciprocal_space_vectors,
    euler_to_matrix,
    sample_to_lab_conversion,
    lab_to_sample_conversion,
)

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
            get_real_space_vectors(a, b, c, alpha, beta, gamma)
        )
        self.a_star_vec_sample, self.b_star_vec_sample, self.c_star_vec_sample = (
            get_reciprocal_space_vectors(a, b, c, alpha, beta, gamma)
        )

        self.a_vec_lab, self.b_vec_lab, self.c_vec_lab = sample_to_lab_conversion(
            self.a_vec_sample, self.b_vec_sample, self.c_vec_sample, roll, pitch, yaw
        )
        self.a_star_vec_lab, self.b_star_vec_lab, self.c_star_vec_lab = (
            sample_to_lab_conversion(
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
