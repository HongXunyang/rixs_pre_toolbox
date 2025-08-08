"""Structure Factor Visualizer for 3D HKL space plotting"""

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np


class StructureFactorVisualizer(FigureCanvas):
    """3D visualizer for structure factors in reciprocal space"""

    def __init__(self, width=8, height=6, dpi=100):
        """Initialize the 3D structure factor visualizer

        Args:
            width: Figure width in inches
            height: Figure height in inches
            dpi: Figure resolution
        """
        # Create figure with proper backend configuration
        self.fig = Figure(figsize=(width, height), dpi=dpi, tight_layout=True)
        super().__init__(self.fig)

        # Initialize state
        self._initialized = False
        self._colorbar = None

        # Create 3D subplot
        self._create_3d_plot()

    def _create_3d_plot(self):
        """Create the 3D subplot and set up basic appearance"""
        try:
            # Clear any existing subplots
            self.fig.clear()

            # Create 3D subplot
            self.axes = self.fig.add_subplot(111, projection="3d")

            # Set basic labels and title
            self.axes.set_xlabel("H (r.l.u.)", fontsize=10)
            self.axes.set_ylabel("K (r.l.u.)", fontsize=10)
            self.axes.set_zlabel("L (r.l.u.)", fontsize=10)
            self.axes.set_title("Structure Factors in Reciprocal Space", fontsize=12)

            # Set default axis limits (0, 2.5)
            self.axes.set_xlim(0, 2.5)
            self.axes.set_ylim(0, 2.5)
            self.axes.set_zlim(0, 2.5)

            # Set integer ticks only
            from matplotlib.ticker import MaxNLocator

            self.axes.xaxis.set_major_locator(MaxNLocator(integer=True))
            self.axes.yaxis.set_major_locator(MaxNLocator(integer=True))
            self.axes.zaxis.set_major_locator(MaxNLocator(integer=True))

            # Set background color
            self.axes.xaxis.pane.fill = False
            self.axes.yaxis.pane.fill = False
            self.axes.zaxis.pane.fill = False

            # Make grid lines less prominent
            self.axes.grid(True, alpha=0.3)

        except Exception as e:
            print(f"Error creating 3D plot: {e}")

    def initialize(self, params: dict = None):
        """Initialize the visualizer

        Args:
            params: Optional parameters (not used for structure factor visualization)
        """
        self._initialized = True
        self._create_3d_plot()
        self.draw()
        return True

    def is_initialized(self):
        """Check if the visualizer is initialized"""
        return self._initialized

    def visualize_structure_factors(self, hkl_list, sf_values):
        """Visualize structure factors as 3D scatter plot

        Args:
            hkl_list: List of [h, k, l] indices, shape (N, 3)
            sf_values: Array of structure factor magnitudes, shape (N,)
        """
        try:
            # Validate inputs
            if len(hkl_list) == 0 or len(sf_values) == 0:
                print("Warning: Empty HKL list or structure factor values")
                return False

            if len(hkl_list) != len(sf_values):
                print(
                    f"Error: HKL list length ({len(hkl_list)}) != SF values length ({len(sf_values)})"
                )
                return False

            # Remove existing colorbar if present
            if self._colorbar is not None:
                self._colorbar.remove()
                self._colorbar = None

            # Clear and recreate the 3D subplot
            self.fig.clear()
            self.axes = self.fig.add_subplot(111, projection="3d")

            # Convert to numpy arrays for easier handling
            hkl_array = np.array(hkl_list)
            sf_array = np.array(sf_values)

            # Extract H, K, L coordinates
            h_coords = hkl_array[:, 0]
            k_coords = hkl_array[:, 1]
            l_coords = hkl_array[:, 2]

            print(f"H coordinates: {h_coords}")
            print(f"K coordinates: {k_coords}")
            print(f"L coordinates: {l_coords}")

            # Calculate dot sizes based on structure factor magnitude
            # Scale factors for visibility
            min_size = 50  # Increased minimum dot size for better visibility
            max_size = 300  # Increased maximum dot size

            print(f"SF array max: {sf_array.max()}, min: {sf_array.min()}")

            # Normalize structure factor values for sizing
            if sf_array.max() > 0:
                normalized_sf = sf_array / sf_array.max()
                dot_sizes = min_size + (max_size - min_size) * normalized_sf
            else:
                dot_sizes = np.full(len(sf_array), min_size)

            print(f"Dot sizes: {dot_sizes}")

            # Create scatter plot
            scatter = self.axes.scatter(
                h_coords,
                k_coords,
                l_coords,
                s=dot_sizes,
                c=sf_array,
                cmap="viridis",
                alpha=0.9,  # More opaque for better visibility
            )

            print(f"Created scatter plot with {len(h_coords)} points")

            # Add text labels next to each dot
            for i, (h, k, l) in enumerate(hkl_list):
                # Create label text (e.g., "010" for [0,1,0])
                label = f"{h}{k}{l}"

                # Add text annotation with slight offset
                self.axes.text(
                    h + 0.1,
                    k + 0.1,
                    l + 0.1,  # Small offset from dot position
                    label,
                    fontsize=8,
                    ha="left",
                    va="bottom",
                    color="black",
                )

            print(f"Added text labels for {len(hkl_list)} points")

            # Add colorbar
            self._colorbar = self.fig.colorbar(
                scatter, ax=self.axes, label="|Structure Factor|", shrink=0.6, pad=0.1
            )

            # Set axis labels and title
            self.axes.set_xlabel("H (r.l.u.)", fontsize=10)
            self.axes.set_ylabel("K (r.l.u.)", fontsize=10)
            self.axes.set_zlabel("L (r.l.u.)", fontsize=10)
            self.axes.set_title("Structure Factors in Reciprocal Space", fontsize=12)

            # Set axis limits with default range (0, 2.5) and adjust if needed
            default_max = 2.5

            # Calculate required ranges based on data
            h_max = max(h_coords.max(), default_max)
            k_max = max(k_coords.max(), default_max)
            l_max = max(l_coords.max(), default_max)

            # Set limits from 0 to calculated max
            self.axes.set_xlim(0, h_max)
            self.axes.set_ylim(0, k_max)
            self.axes.set_zlim(0, l_max)

            # Set integer ticks only
            from matplotlib.ticker import MaxNLocator

            self.axes.xaxis.set_major_locator(MaxNLocator(integer=True))
            self.axes.yaxis.set_major_locator(MaxNLocator(integer=True))
            self.axes.zaxis.set_major_locator(MaxNLocator(integer=True))

            print(f"Set axis limits: X(0, {h_max}), Y(0, {k_max}), Z(0, {l_max})")

            # Improve 3D viewing angle
            self.axes.view_init(elev=20, azim=45)

            # Enable grid with low alpha
            self.axes.grid(True, alpha=0.3)

            # Update the canvas
            self.draw()

            print(f"Successfully plotted {len(hkl_list)} structure factors")
            return True

        except Exception as e:
            print(f"Error in visualize_structure_factors: {e}")
            return False

    def clear_plot(self):
        """Clear the current plot"""
        try:
            # Remove existing colorbar if present
            if self._colorbar is not None:
                self._colorbar.remove()
                self._colorbar = None

            # Recreate the 3D plot
            self._create_3d_plot()
            self.draw()
        except Exception as e:
            print(f"Error clearing plot: {e}")

    def visualize(self):
        """Create an empty visualization (placeholder)"""
        try:
            self._create_3d_plot()
            self.draw()
            return True
        except Exception as e:
            print(f"Error in visualize: {e}")
            return False
