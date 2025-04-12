#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from matplotlib.figure import Figure


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

        # Draw simple cube as placeholder for Brillouin zone
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
                alpha=0.3,
            )

        # Draw origin
        ax.scatter([0], [0], [0], color="red", s=100)

        # Draw reciprocal axes
        ax.quiver(0, 0, 0, 1, 0, 0, color="r", arrow_length_ratio=0.1, label="a*")
        ax.quiver(0, 0, 0, 0, 1, 0, color="g", arrow_length_ratio=0.1, label="b*")
        ax.quiver(0, 0, 0, 0, 0, 1, color="b", arrow_length_ratio=0.1, label="c*")

        # Draw Q vector if present in data
        if "h" in data and "k" in data and "l" in data:
            h, k, l = data["h"], data["k"], data["l"]

            # Scale for visibility
            scale = 1.5 / max(1e-6, np.sqrt(h**2 + k**2 + l**2))

            # Draw Q vector
            ax.quiver(
                0,
                0,
                0,
                h * scale,
                k * scale,
                l * scale,
                color="m",
                arrow_length_ratio=0.1,
                linewidth=2,
            )

            # Draw HKL point
            ax.scatter([h * scale], [k * scale], [l * scale], color="purple", s=100)

            # Annotate the point
            ax.text(
                h * scale,
                k * scale,
                l * scale,
                f"({h:.2f}, {k:.2f}, {l:.2f})",
                color="purple",
            )

        # Set labels and title
        ax.set_xlabel("a* (Å⁻¹)")
        ax.set_ylabel("b* (Å⁻¹)")
        ax.set_zlabel("c* (Å⁻¹)")
        ax.set_title("Scattering Geometry")

        # Set equal aspect ratio
        ax.set_box_aspect([1, 1, 1])

        return fig
