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
from packages.classes import Lattice

class Sample:
    """This is a class for the sample."""

    def __init__(self):
        """Initialize the sample."""
        self.lattice = Lattice()
        self.roll = 0
        self.pitch = 0
        self.yaw = 0
        self.a_vec_sample = None
        self.b_vec_sample = None
        self.c_vec_sample = None
        self.a_star_vec_sample = None
        self.b_star_vec_sample = None
        self.c_star_vec_sample = None

    def initialize(self, a, b, c, alpha, beta, gamma, roll, pitch, yaw):
        """Initialize the sample.

        Args:
            a, b, c (float): Lattice constants in Angstroms
            alpha, beta, gamma (float): Lattice angles in degrees
            roll, pitch, yaw (float): Euler angles in degrees
        """
        self.lattice.initialize(a, b, c, alpha, beta, gamma)
        self.roll = roll
        self.pitch = pitch
        self.yaw = yaw
        self.calculate_real_space_vectors()
        self.calculate_reciprocal_space_vectors()

    def get_lattice_parameters(self):
        """Get the parameters of the sample. if None, raise "initialize the sample first"."""
        try:
            return self.lattice.get_lattice_parameters()
        except KeyError as exc:
            raise ValueError("initialize the sample first") from exc

    def get_real_space_vectors(self):
        """Get the real space vectors in the sample frame."""
        return self.a_vec_sample, self.b_vec_sample, self.c_vec_sample

    def get_reciprocal_space_vectors(self):
        """Get the reciprocal space vectors in the sample frame."""
        return self.a_star_vec_sample, self.b_star_vec_sample, self.c_star_vec_sample

    def calculate_real_space_vectors(self):
        """Get the real space vectors in the sample frame."""
        a_vec_lattice, b_vec_lattice, c_vec_lattice = (
            self.lattice.get_real_space_vectors()
        )
        rotation_matrix = euler_to_matrix(self.roll, self.pitch, self.yaw)
        self.a_vec_sample = rotation_matrix @ a_vec_lattice
        self.b_vec_sample = rotation_matrix @ b_vec_lattice
        self.c_vec_sample = rotation_matrix @ c_vec_lattice

    def calculate_reciprocal_space_vectors(self):
        """Get the reciprocal space vectors in the sample frame."""
        a_star_vec_lattice, b_star_vec_lattice, c_star_vec_lattice = (
            self.lattice.get_reciprocal_space_vectors()
        )
        rotation_matrix = euler_to_matrix(self.roll, self.pitch, self.yaw)
        self.a_star_vec_sample = rotation_matrix @ a_star_vec_lattice
        self.b_star_vec_sample = rotation_matrix @ b_star_vec_lattice
        self.c_star_vec_sample = rotation_matrix @ c_star_vec_lattice
