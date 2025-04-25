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
        self.axes.view_init(elev=25, azim=-67)

        # Remove ticks and grid
        self.axes.set_xticks([])
        self.axes.set_yticks([])
        self.axes.set_zticks([])
        self.axes.grid(False)

        # Set initial limits
        self.axes.set_xlim(-0.75, 0.75)
        self.axes.set_ylim(-0.75, 0.75)
        self.axes.set_zlim(-0.75, 0.75)

    def visualize_scattering_geometry(self, scattering_angles):
        """Update the visualization with new scattering angles."""
        # Clear previous plot
        self.axes.clear()

        # Define origin and cube vertices
        origin = np.array([-0.5, -0.125, -0.5])

        # Define vertices of a unit cube
        ver = np.array(
            [
                [1, 1, 0],
                [0, 1, 0],
                [0, 1, 1],
                [1, 1, 1],
                [0, 0, 1],
                [1, 0, 1],
                [1, 0, 0],
                [0, 0, 0],
            ]
        )

        # Define faces of the unit cube
        fac = np.array(
            [
                [1, 2, 3, 4],
                [4, 3, 5, 6],
                [6, 7, 8, 5],
                [1, 2, 8, 7],
                [6, 7, 1, 4],
                [2, 3, 5, 8],
            ]
        )

        # Cube dimensions and color
        X, Y, Z = 1, 0.25, 1
        color = [0.3, 0.3, 0.3]

        # Scale vertices by dimensions and translate to origin
        cube = np.array(
            [
                ver[:, 0] * X + origin[0],
                ver[:, 1] * Y + origin[1] + Y / 2,
                ver[:, 2] * Z + origin[2],
            ]
        ).T

        # Plot the cube
        self.axes.add_collection3d(
            Poly3DCollection(
                cube[fac - 1], facecolors=color, edgecolors="black", alpha=0.3
            )
        )

        # Plot scattering plane
        vertices = np.array(
            [
                [origin[0] - 0.25, origin[1] + 0.375, 0],
                [origin[0] + 1.25, origin[1] + 0.375, 0],
                [origin[0] + 1.25, origin[1] - 0.75, 0],
                [origin[0] - 0.25, origin[1] - 0.75, 0],
            ]
        )

        faces = np.array([[0, 1, 2, 3]])
        self.axes.add_collection3d(
            Poly3DCollection(
                vertices[faces],
                facecolors=[0.3010, 0.7450, 0.9330],
                edgecolors=[0.7, 0.7, 0.7],
                alpha=0.15,
            )
        )

        # Extract angles from data
        theta = scattering_angles.get("theta", 0)  # theta angle
        tth = scattering_angles.get("tth", 0)  # two theta angle

        # Plot incident beam (k_in)
        k_in_length = 0.8
        k_in_x = k_in_length * np.cos(np.radians(theta))
        k_in_y = -k_in_length * np.sin(np.radians(theta))
        k_in_z = 0
        # Draw colored arrow on top
        self.axes.quiver(
            k_in_x,
            k_in_y,
            0,
            -k_in_x,
            -k_in_y,
            0,
            color=[0.93, 0.694, 0.125],
            alpha=0.8,
            linewidth=4,
            arrow_length_ratio=0.3,
            zorder=10,
        )

        # Plot scattered beam (k_out)
        k_out_length = 0.8
        k_out_x = -k_out_length * np.cos(np.radians(tth - theta))
        k_out_y = -k_out_length * np.sin(np.radians(tth - theta))
        k_out_z = 0

        # Draw colored arrow on top
        self.axes.quiver(
            0,
            0,
            0,
            k_out_x,
            k_out_y,
            k_out_z,
            color=[0, 0.5, 0],
            alpha=0.8,
            linewidth=4,
            arrow_length_ratio=0.3,
            zorder=10,
        )

        # remove ticks
        self.axes.set_xticks([])
        self.axes.set_yticks([])
        self.axes.set_zticks([])

        # Update the canvas
        self.draw()
