#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


class BrillouinVisualizer:
    """Class for visualizing Brillouin zone and scattering geometry."""

    def __init__(self, reciprocal_lattice):
        """Initialize the visualizer with reciprocal lattice parameters.

        Args:
            reciprocal_lattice (dict): Dictionary containing reciprocal lattice vectors
        """
        self.reciprocal_lattice = reciprocal_lattice

    def visualize_brillouin_zone(self):
        """Visualize the Brillouin zone.

        Returns:
            Figure: Matplotlib figure with the visualization
        """
        # Create a figure for the Brillouin zone
        fig = Figure(figsize=(6, 5))
        ax = fig.add_subplot(111, projection="3d")

        # Draw simple cube as placeholder for Brillouin zone
        # In a real implementation, this would calculate and draw
        # the actual Brillouin zone for the crystal structure
        r = 1.0

        # Vertices of a cube
        vertices = np.array(
            [
                [r, r, r],
                [r, r, -r],
                [r, -r, r],
                [r, -r, -r],
                [-r, r, r],
                [-r, r, -r],
                [-r, -r, r],
                [-r, -r, -r],
            ]
        )

        # Draw edges
        edges = [
            (0, 1),
            (0, 2),
            (0, 4),
            (1, 3),
            (1, 5),
            (2, 3),
            (2, 6),
            (3, 7),
            (4, 5),
            (4, 6),
            (5, 7),
            (6, 7),
        ]

        for edge in edges:
            ax.plot3D(
                [vertices[edge[0], 0], vertices[edge[1], 0]],
                [vertices[edge[0], 1], vertices[edge[1], 1]],
                [vertices[edge[0], 2], vertices[edge[1], 2]],
                "b-",
            )

        # Draw origin
        ax.scatter([0], [0], [0], color="red", s=100)

        # Draw reciprocal axes
        ax.quiver(0, 0, 0, 1, 0, 0, color="r", arrow_length_ratio=0.1, label="a*")
        ax.quiver(0, 0, 0, 0, 1, 0, color="g", arrow_length_ratio=0.1, label="b*")
        ax.quiver(0, 0, 0, 0, 0, 1, color="b", arrow_length_ratio=0.1, label="c*")

        # Set labels and title
        ax.set_xlabel("a* (Å⁻¹)")
        ax.set_ylabel("b* (Å⁻¹)")
        ax.set_zlabel("c* (Å⁻¹)")
        ax.set_title("Brillouin Zone")

        # Set equal aspect ratio
        ax.set_box_aspect([1, 1, 1])

        return fig

    def visualize_scattering_geometry(self, data):
        """Visualize the scattering geometry.

        Args:
            data (dict): Data containing HKL indices and scattering angles

        Returns:
            Figure: Matplotlib figure with the visualization
        """
        # Create a figure for the scattering geometry
        fig = Figure(figsize=(6, 5))
        ax = fig.add_subplot(111, projection="3d")

        # Set 3D view angle
        ax.view_init(elev=25, azim=-67)

        # Remove ticks and grid
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])
        ax.grid(False)

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
                ver[:, 1] * Y + origin[1],
                ver[:, 2] * Z + origin[2],
            ]
        ).T

        # Plot the cube
        ax.add_collection3d(
            Poly3DCollection(
                cube[fac - 1], facecolors=color, edgecolors="black", alpha=0.8
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
        ax.add_collection3d(
            Poly3DCollection(
                vertices[faces],
                facecolors=[0.3010, 0.7450, 0.9330],
                edgecolors=[0.7, 0.7, 0.7],
                alpha=0.3,
            )
        )

        # Extract angles from data
        th = data.get("th", 0)  # theta angle
        tth = data.get("tth", 0)  # two theta angle

        # Plot incident beam (k_in)
        k_in_length = 0.8
        k_in_x = k_in_length * np.cos(np.radians(th))
        k_in_y = -k_in_length * np.sin(np.radians(th))
        k_in_z = 0

        # Draw white outline first
        ax.quiver(
            k_in_x,
            k_in_y,
            0,
            -k_in_x,
            -k_in_y,
            0,
            color="white",
            alpha=0.8,
            linewidth=4,
            arrow_length_ratio=0.3,
        )
        # Draw colored arrow on top
        ax.quiver(
            k_in_x,
            k_in_y,
            0,
            -k_in_x,
            -k_in_y,
            0,
            color=[0.93, 0.694, 0.125],
            alpha=0.8,
            linewidth=2,
            arrow_length_ratio=0.3,
        )

        # Plot scattered beam (k_out)
        k_out_length = 0.8
        k_out_x = -k_out_length * np.cos(np.radians(tth - th))
        k_out_y = -k_out_length * np.sin(np.radians(tth - th)) - 0.125
        k_out_z = 0

        # Draw white outline first
        ax.quiver(
            0,
            -0.125,
            0,
            k_out_x,
            k_out_y,
            k_out_z,
            color="white",
            alpha=0.8,
            linewidth=4,
            arrow_length_ratio=0.3,
        )
        # Draw colored arrow on top
        ax.quiver(
            0,
            -0.125,
            0,
            k_out_x,
            k_out_y,
            k_out_z,
            color=[0, 0.5, 0],
            alpha=0.8,
            linewidth=2,
            arrow_length_ratio=0.3,
        )

        # Set axis properties - zoomed in view
        ax.set_xlim(-0.75, 0.75)
        ax.set_ylim(-0.75, 0.75)
        ax.set_zlim(-0.75, 0.75)

        # Add title with HKL values if present
        h = data.get("h", 0)
        k = data.get("k", 0)
        l = data.get("l", 0)
        ax.set_title(f"H={h:.3f}, K={k:.3f}, L={l:.3f}", pad=20, fontsize=14)

        return fig
