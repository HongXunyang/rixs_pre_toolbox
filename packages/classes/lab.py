"""This is a class for the lab."""

import sys
import os
import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from packages.classes import Sample
from packages.utils import angle_to_matrix


class Lab:
    """This is a class for the lab."""

    def __init__(self):
        """Initialize the lab."""
        self.sample = Sample()
        self.theta = 0
        self.phi = 0
        self.chi = 0
        self.a_vec_lab = None
        self.b_vec_lab = None
        self.c_vec_lab = None
        self.a_star_vec_lab = None
        self.b_star_vec_lab = None
        self.c_star_vec_lab = None

    def initialize(
        self, a, b, c, alpha, beta, gamma, roll, pitch, yaw, theta, phi, chi
    ):
        """Initialize the lab."""
        self.sample.initialize(a, b, c, alpha, beta, gamma, roll, pitch, yaw)
        self.theta = theta
        self.phi = phi
        self.chi = chi
        self.calculate_real_space_vectors()
        self.calculate_reciprocal_space_vectors()

    def get_lattice_parameters(self):
        """Get the parameters of the sample."""
        return self.sample.get_lattice_parameters()

    def get_real_space_vectors(self):
        """Get the real space vectors in the lab frame."""
        return self.a_vec_lab, self.b_vec_lab, self.c_vec_lab

    def get_reciprocal_space_vectors(self):
        """Get the reciprocal space vectors in the lab frame."""
        return self.a_star_vec_lab, self.b_star_vec_lab, self.c_star_vec_lab

    def calculate_real_space_vectors(self):
        """Get the real space vectors in the lab frame."""
        a_vec_sample, b_vec_sample, c_vec_sample = self.sample.get_real_space_vectors()
        rotation_matrix = angle_to_matrix(self.theta, self.phi, self.chi)
        self.a_vec_lab = rotation_matrix @ a_vec_sample
        self.b_vec_lab = rotation_matrix @ b_vec_sample
        self.c_vec_lab = rotation_matrix @ c_vec_sample

    def calculate_reciprocal_space_vectors(self):
        """Get the reciprocal space vectors in the lab frame."""
        a_star_vec_sample, b_star_vec_sample, c_star_vec_sample = (
            self.sample.get_reciprocal_space_vectors()
        )
        rotation_matrix = angle_to_matrix(self.theta, self.phi, self.chi)
        self.a_star_vec_lab = rotation_matrix @ a_star_vec_sample
        self.b_star_vec_lab = rotation_matrix @ b_star_vec_sample
        self.c_star_vec_lab = rotation_matrix @ c_star_vec_sample
