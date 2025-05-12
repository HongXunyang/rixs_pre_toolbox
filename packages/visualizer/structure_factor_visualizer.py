"""This is a class to visualize the structure factor"""

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np


class StructureFactorVisualizer(FigureCanvas):
    """Visualizer for structure factor"""

    def __init__(self, width=6, height=6, dpi=100):
        """Initialize the visualizer"""
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111, projection="3d")
        super().__init__(self.fig)
        self._initialized = False

    def initialize(self, params: dict):
        """Initialize the visualizer with lattice parameters"""
        self._initialized = True
        return True

    def is_initialized(self):
        """Check if the visualizer is initialized"""
        return self._initialized

    def visualize_structure_factors(self, hkl_list, sf_values):
        """Visualize the structure factors in 3D space

        Args:
            hkl_list: List of (h, k, l) tuples
            sf_values: List of structure factor values
        """
        # Clear the figure
        self.axes.clear()

        # Extract h, k, l values for plotting
        h_values = [hkl[0] for hkl in hkl_list]
        k_values = [hkl[1] for hkl in hkl_list]
        l_values = [hkl[2] for hkl in hkl_list]

        # Plot the structure factors
        scatter = self.axes.scatter(
            h_values,
            k_values,
            l_values,
            s=[v * 50 for v in sf_values],  # Size based on structure factor
            c=sf_values,  # Color based on structure factor
            cmap="viridis",
            alpha=0.6,
        )

        # Add a colorbar
        self.fig.colorbar(scatter, ax=self.axes, label="Structure Factor")

        # Set labels and limits
        self.axes.set_xlabel("H (r.l.u.)")
        self.axes.set_ylabel("K (r.l.u.)")
        self.axes.set_zlabel("L (r.l.u.)")
        self.axes.set_title("Structure Factors in Reciprocal Space")

        # Draw the updated figure
        self.draw()

        return True

    def visualize(self):
        """Visualize the structure factor (placeholder)"""
        self.axes.clear()
        self.axes.set_xlabel("H (r.l.u.)")
        self.axes.set_ylabel("K (r.l.u.)")
        self.axes.set_zlabel("L (r.l.u.)")
        self.axes.set_title("Structure Factors in Reciprocal Space")
        self.draw()
        return True

    def update(self):
        """Update the visualizer (placeholder)"""
        self.draw()
        return True
