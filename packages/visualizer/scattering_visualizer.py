#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""This is a class to visualize the X-ray scattering geometry"""

import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


class ScatteringVisualizer(FigureCanvas):
    """Visualizer for scattering geometry with 3D interactive canvas."""

    def __init__(self, width=4, height=4, dpi=100):
        """Initialize the visualizer with a 3D canvas."""
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111, projection="3d")
        super().__init__(self.fig)

        # Set initial view
        self.axes.view_init(elev=39, azim=75)

        self.e_H = np.array([1, 0, 0])
        self.e_K = np.array([0, 1, 0])
        self.e_L = np.array([0, 0, 1])

    def initialize(self, e_H, e_K, e_L):
        """Initialize the visualizer with the given crystal coordinates system."""
        self.e_H = e_H
        self.e_K = e_K
        self.e_L = e_L
        return True

    def visualize_lab_system(self, chi=0, phi=0, is_clear=True):
        """Update the visualization with new crystal coordinates system."""
        if is_clear:
            # Clear previous plot
            self.axes.clear()

        # Define vertices of the sample
        vertices_sample = np.array(
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
        vertices_sample = _rotate_vertices(vertices_sample, phi, chi)
        # Define faces of the cube
        faces_sample = np.array(
            [
                [0, 1, 2, 3],  # top face
                [4, 5, 6, 7],  # bottom face
                [0, 1, 5, 4],  # front face
                [2, 3, 7, 6],  # back face
                [0, 3, 7, 4],  # right face
                [1, 2, 6, 5],  # left face
            ]
        )

        # Plot the cube
        self.axes.add_collection3d(
            Poly3DCollection(
                vertices_sample[faces_sample],
                facecolors=[0.3, 0.3, 0.3],
                edgecolors=[0.55, 0.55, 0.55],
                alpha=0.2,
            )
        )

        # Normalize the vectors
        e_H_norm = self.e_H / np.linalg.norm(self.e_H)
        e_K_norm = self.e_K / np.linalg.norm(self.e_K)
        e_L_norm = self.e_L / np.linalg.norm(self.e_L)

        # Plot the normalized vectors
        vectors = [e_H_norm, e_K_norm, e_L_norm]
        vectors = _rotate_vertices(vectors, phi, chi)
        colors = ["r", "g", "b"]
        labels = ["$e_H$", "$e_K$", "$e_L$"]

        for vec, color, label in zip(vectors, colors, labels):
            # Plot the vector
            self.axes.quiver(
                0,
                0,
                0,  # origin
                vec[0],
                vec[1],
                vec[2],  # vector components
                color=color,
                alpha=0.25,
                linewidth=2,
                arrow_length_ratio=0.2,
            )

            # Add text label at the tip of the vector
            # Add a small offset to prevent text from overlapping with the arrow
            offset = 0.3
            self.axes.text(
                vec[0] + offset,
                vec[1] + offset,
                vec[2] + offset,
                label,
                color=color,
                fontsize=14,
                ha="center",
                alpha=0.25,
            )

        # plot the vector of the lab coordinate system, by default it is the unit vectors
        e_X = np.array([1, 0, 0]) / 0.65
        e_Y = np.array([0, 1, 0]) / 0.65
        e_Z = np.array([0, 0, 1]) / 0.65
        vectors = [e_X, e_Y, e_Z]
        vectors = _rotate_vertices(vectors, phi, chi)
        colors = ["r", "g", "b"]
        labels = ["$X$", "$Y$", "$Z$"]

        for vec, color, label in zip(vectors, colors, labels):
            # Plot the vector
            self.axes.quiver(
                0,
                0,
                0,  # origin
                vec[0],
                vec[1],
                vec[2],  # vector components
                color=(64 / 255, 148 / 255, 184 / 255),
                alpha=0.4,
                linewidth=0.8,
                arrow_length_ratio=0.1,
            )

            # Add text label at the tip of the vector
            # Add a small offset to prevent text from overlapping with the arrow
            offset = 0.2
            self.axes.text(
                vec[0] + offset,
                vec[1] + offset,
                vec[2] + offset,
                label,
                color=(64 / 255, 148 / 255, 184 / 255),
                alpha=0.4,
            )

        # Update the canvas
        # Set axis limits
        self.axes.set_xlim(-1, 1)
        self.axes.set_ylim(-1, 1)
        self.axes.set_zlim(-1, 1)

        # Remove ticks
        self.axes.set_xticks([])
        self.axes.set_yticks([])
        self.axes.set_zticks([])

        self.draw()

    def visualize_scattering_geometry(self, scattering_angles=None, is_clear=True):
        """Update the visualization with new scattering angles."""
        if is_clear:
            # Clear previous plot
            self.axes.clear()

        # Plot the scattering plane
        scatter_plane_vertices = np.array(
            [
                [0.75, 0, -0.25],  # bottom right
                [-0.75, 0, -0.25],  # bottom left
                [0.75, 0, 1.25],  # top right
                [-0.75, 0, 1.25],  # top left
            ]
        )

        scatter_plane_faces = np.array([[0, 1, 3, 2]])  # single face
        self.axes.add_collection3d(
            Poly3DCollection(
                scatter_plane_vertices[scatter_plane_faces],
                facecolors=[0.3010, 0.7450, 0.9330],  # light blue
                edgecolors=[0.7, 0.7, 0.7],
                alpha=0.15,
            )
        )

        # Plot the x-ray beam
        if scattering_angles is None:
            scattering_angles = {
                "theta": 50,
                "tth": 150,
            }

        # Extract angles from data
        theta = scattering_angles.get("theta", 50)  # theta angle
        tth = scattering_angles.get("tth", 150)  # two theta angle

        # Plot incident beam (k_in)
        offset = 0
        k_in_length = 1.3
        k_in_x = k_in_length * np.cos(np.radians(theta))
        k_in_z = k_in_length * np.sin(np.radians(theta))
        k_in_y = 0
        # Draw colored arrow on top
        self.axes.quiver(
            k_in_x,
            k_in_y + offset,
            k_in_z,
            -k_in_x,
            -k_in_y,
            -k_in_z,
            color=(191 / 255, 44 / 255, 0),
            alpha=1,
            linewidth=5,
            arrow_length_ratio=0.2,
            zorder=10,
        )

        # Plot scattered beam (k_out)
        k_out_length = 1.3
        k_out_x = -k_out_length * np.cos(np.radians(tth - theta))
        k_out_z = k_out_length * np.sin(np.radians(tth - theta))
        k_out_y = 0

        # Draw colored arrow on top
        self.axes.quiver(
            0,
            0 + offset,
            0,
            k_out_x,
            k_out_y,
            k_out_z,
            color=(2 / 255, 78 / 255, 191 / 255),
            linewidth=5,
            arrow_length_ratio=0.2,
            zorder=10,
        )

        # Set axis limits
        self.axes.set_xlim(-1, 1)
        self.axes.set_ylim(-1, 1)
        self.axes.set_zlim(-1, 1)

        # Remove ticks
        self.axes.set_xticks([])
        self.axes.set_yticks([])
        self.axes.set_zticks([])
        # Update the canvas
        self.draw()


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
