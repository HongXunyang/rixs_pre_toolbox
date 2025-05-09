"""This is the sample class for the sample object."""

import numpy as np
from copy import deepcopy


class Sample:
    """This is the sample class for the sample object."""

    def __init__(self, parent_tab=None):
        """Initialize the sample object."""
        self.e_H = np.array([1, 0, 0])
        self.e_K = np.array([0, 1, 0])
        self.e_L = np.array([0, 0, 1])
        self.e_X = np.array([1, 0, 0])
        self.e_Y = np.array([0, 1, 0])
        self.e_Z = np.array([0, 0, 1])
        self.unit_vectors_momentum = [self.e_H, self.e_K, self.e_L]
        self.unit_vectors_coordinates = [self.e_X, self.e_Y, self.e_Z]
        self.theta = 0
        self.phi = 0
        self.chi = 0
        # Define default vertices of the sample
        self.vertices = np.array(
            [
                [0.5, 0.5, 0],  # top front right
                [0.5, -0.5, 0],  # top front left
                [-0.5, -0.5, 0],  # top back left
                [-0.5, 0.5, 0],  # top back right
                [0.5, 0.5, -0.25],  # bottom front right
                [0.5, -0.5, -0.25],  # bottom front left
                [-0.5, -0.5, -0.25],  # bottom back left
                [-0.5, 0.5, -0.25],  # bottom back right
            ]
        )
        self.vertices_default = deepcopy(self.vertices)
        # Define default faces of the cube
        self.faces = np.array(
            [
                [0, 1, 2, 3],  # top face
                [4, 5, 6, 7],  # bottom face
                [0, 1, 5, 4],  # front face
                [2, 3, 7, 6],  # back face
                [0, 3, 7, 4],  # right face
                [1, 2, 6, 5],  # left face
            ]
        )

    def initialize(self, e_H, e_K, e_L):
        """Initialize the sample object."""
        self.e_H = e_H / np.linalg.norm(e_H)
        self.e_K = e_K / np.linalg.norm(e_K)
        self.e_L = e_L / np.linalg.norm(e_L)
        self.unit_vectors_momentum = [self.e_H, self.e_K, self.e_L]

    def rotate_sample(self, theta, phi, chi):
        """Rotate the sample object."""
        self.theta = theta
        self.phi = phi
        self.chi = chi
        self.vertices = _rotate_vertices(self.vertices_default, phi, chi)
        self.unit_vectors_coordinates = _rotate_vertices(
            self.unit_vectors_coordinates, phi, chi
        )
        self.unit_vectors_momentum = _rotate_vertices(
            self.unit_vectors_momentum, phi, chi
        )
        self.e_H, self.e_K, self.e_L = self.unit_vectors_momentum
        self.e_X, self.e_Y, self.e_Z = self.unit_vectors_coordinates

    def get_rearranged_vertices(self):
        """Get the rearranged vertices of the sample for 3D plotting."""
        return self.vertices[self.faces]

    def get_sample_orientation(self):
        """Get the sample orientation."""
        return self.theta, self.phi, self.chi


def _rotate_vertices(vertices, phi, chi):
    """Rotate the vertices of the cube by the given phi and chi angles."""
    # Convert angles to radians
    phi_rad = np.radians(phi)
    chi_rad = np.radians(chi)

    chi_mat_sample = np.array(
        [
            [1, 0, 0],
            [0, np.cos(chi_rad), -np.sin(chi_rad)],
            [0, np.sin(chi_rad), np.cos(chi_rad)],
        ]
    )
    phi_mat_sample = np.array(
        [
            [np.cos(phi_rad), -np.sin(phi_rad), 0],
            [np.sin(phi_rad), np.cos(phi_rad), 0],
            [0, 0, 1],
        ]
    )
    vertices = np.array(vertices)
    for i, vertex in enumerate(vertices):
        vertex_temp = chi_mat_sample @ phi_mat_sample @ vertex
        vertices[i] = vertex_temp
    return vertices
