#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from matplotlib.figure import Figure


class BrillouinVisualizer:
    """Interface for the Brillouin zone visualizer.
    
    This class handles all the visualization functionality for the Brillouin zone
    visualizer tab. It's a pure Python implementation without PyQt dependencies.
    """
    
    def __init__(self):
        """Initialize the visualizer."""
        self._initialized = False
        self.crystal_structure = None
        self.brillouin_zone = None
        self.special_points = {}
    
    def initialize_from_cif(self, cif_file_path):
        """Initialize the visualizer from a CIF file.
        
        Args:
            cif_file_path (str): Path to the CIF file
            
        Returns:
            bool: True if initialization was successful
        """
        try:
            # Placeholder implementation
            # In a real implementation, this would parse the CIF file
            # and calculate the Brillouin zone and special points
            
            # Mock data for demonstration
            self.crystal_structure = {
                'name': 'Silicon',
                'space_group': 'Fd-3m',
                'lattice_type': 'FCC',
                'lattice_constants': {
                    'a': 5.43,
                    'b': 5.43,
                    'c': 5.43,
                    'alpha': 90.0,
                    'beta': 90.0,
                    'gamma': 90.0
                }
            }
            
            # Mock Brillouin zone (just vertices for a cube)
            self.brillouin_zone = {
                'vertices': np.array([
                    [1, 1, 1], [1, 1, -1], [1, -1, 1], [1, -1, -1],
                    [-1, 1, 1], [-1, 1, -1], [-1, -1, 1], [-1, -1, -1]
                ]),
                'edges': [
                    (0, 1), (0, 2), (0, 4), (1, 3), (1, 5), (2, 3),
                    (2, 6), (3, 7), (4, 5), (4, 6), (5, 7), (6, 7)
                ]
            }
            
            # Mock special points for FCC Brillouin zone
            self.special_points = {
                'Γ': [0, 0, 0],
                'X': [0, 0, 1],
                'L': [0.5, 0.5, 0.5],
                'W': [0.5, 0.25, 0.75],
                'K': [0.375, 0.375, 0.75]
            }
            
            self._initialized = True
            return True
        except Exception as e:
            print(f"Error initializing visualizer from CIF file: {str(e)}")
            return False
    
    def initialize_from_parameters(self, lattice_type, lattice_constants):
        """Initialize visualizer from lattice parameters.
        
        Args:
            lattice_type (str): Type of lattice ('SC', 'BCC', 'FCC', etc.)
            lattice_constants (dict): Dictionary of lattice constants
            
        Returns:
            bool: True if initialization was successful
        """
        try:
            # Placeholder implementation
            # In a real implementation, this would calculate the Brillouin zone
            # and special points from the lattice parameters
            
            # Mock data
            self.crystal_structure = {
                'name': 'Custom',
                'space_group': 'P1',
                'lattice_type': lattice_type,
                'lattice_constants': lattice_constants
            }
            
            # For simplicity, using same mock Brillouin zone as above
            # In a real implementation, this would be calculated properly
            self.brillouin_zone = {
                'vertices': np.array([
                    [1, 1, 1], [1, 1, -1], [1, -1, 1], [1, -1, -1],
                    [-1, 1, 1], [-1, 1, -1], [-1, -1, 1], [-1, -1, -1]
                ]),
                'edges': [
                    (0, 1), (0, 2), (0, 4), (1, 3), (1, 5), (2, 3),
                    (2, 6), (3, 7), (4, 5), (4, 6), (5, 7), (6, 7)
                ]
            }
            
            # Set special points based on lattice type
            if lattice_type == 'SC':
                self.special_points = {
                    'Γ': [0, 0, 0],
                    'X': [0.5, 0, 0],
                    'M': [0.5, 0.5, 0],
                    'R': [0.5, 0.5, 0.5]
                }
            elif lattice_type == 'BCC':
                self.special_points = {
                    'Γ': [0, 0, 0],
                    'H': [0.5, -0.5, 0.5],
                    'N': [0, 0, 0.5],
                    'P': [0.25, 0.25, 0.25]
                }
            elif lattice_type == 'FCC':
                self.special_points = {
                    'Γ': [0, 0, 0],
                    'X': [0, 0, 1],
                    'L': [0.5, 0.5, 0.5],
                    'W': [0.5, 0.25, 0.75],
                    'K': [0.375, 0.375, 0.75]
                }
            else:
                # Default simple points
                self.special_points = {
                    'Γ': [0, 0, 0],
                    'X': [0.5, 0, 0],
                    'Y': [0, 0.5, 0],
                    'Z': [0, 0, 0.5]
                }
            
            self._initialized = True
            return True
        except Exception as e:
            print(f"Error initializing visualizer from parameters: {str(e)}")
            return False
    
    def is_initialized(self):
        """Check if the visualizer is initialized.
        
        Returns:
            bool: True if the visualizer is initialized
        """
        return self._initialized
    
    def get_crystal_structure(self):
        """Get the crystal structure information.
        
        Returns:
            dict: Dictionary containing crystal structure information
        """
        return self.crystal_structure
    
    def get_special_points(self):
        """Get the special points in the Brillouin zone.
        
        Returns:
            dict: Dictionary mapping point names to coordinates
        """
        return self.special_points
    
    def visualize_brillouin_zone(self, show_special_points=True, view_direction=None):
        """Visualize the Brillouin zone.
        
        Args:
            show_special_points (bool): Whether to show special points
            view_direction (str): View direction ('x', 'y', 'z' or None for 3D)
            
        Returns:
            Figure: Matplotlib figure with the visualization
        """
        if not self.is_initialized():
            return None
        
        # Create figure
        fig = Figure(figsize=(6, 5))
        
        # Create 3D axes or 2D projection based on view_direction
        if view_direction is None:
            ax = fig.add_subplot(111, projection='3d')
            is_3d = True
        else:
            ax = fig.add_subplot(111)
            is_3d = False
        
        # Draw Brillouin zone
        vertices = self.brillouin_zone['vertices']
        edges = self.brillouin_zone['edges']
        
        if is_3d:
            # Draw edges in 3D
            for edge in edges:
                ax.plot3D(
                    [vertices[edge[0], 0], vertices[edge[1], 0]],
                    [vertices[edge[0], 1], vertices[edge[1], 1]],
                    [vertices[edge[0], 2], vertices[edge[1], 2]],
                    'b-'
                )
            
            # Draw origin
            ax.scatter([0], [0], [0], color='red', s=100)
            
            # Draw special points if requested
            if show_special_points and self.special_points:
                for name, coords in self.special_points.items():
                    ax.scatter(coords[0], coords[1], coords[2], color='green', s=50)
                    ax.text(coords[0], coords[1], coords[2], name, fontsize=10)
            
            # Set labels
            ax.set_xlabel('kx')
            ax.set_ylabel('ky')
            ax.set_zlabel('kz')
            
            # Set equal aspect ratio
            ax.set_box_aspect([1, 1, 1])
        else:
            # Create 2D projection
            if view_direction == 'x':
                idx1, idx2 = 1, 2  # y, z
                xlabel, ylabel = 'ky', 'kz'
            elif view_direction == 'y':
                idx1, idx2 = 0, 2  # x, z
                xlabel, ylabel = 'kx', 'kz'
            else:  # 'z' or default
                idx1, idx2 = 0, 1  # x, y
                xlabel, ylabel = 'kx', 'ky'
            
            # Draw projection of edges
            for edge in edges:
                ax.plot(
                    [vertices[edge[0], idx1], vertices[edge[1], idx1]],
                    [vertices[edge[0], idx2], vertices[edge[1], idx2]],
                    'b-'
                )
            
            # Draw origin
            ax.scatter([0], [0], color='red', s=100)
            
            # Draw special points if requested
            if show_special_points and self.special_points:
                for name, coords in self.special_points.items():
                    ax.scatter(coords[idx1], coords[idx2], color='green', s=50)
                    ax.text(coords[idx1], coords[idx2], name, fontsize=10)
            
            # Set labels
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            
            # Set equal aspect ratio
            ax.set_aspect('equal')
        
        # Set title
        if self.crystal_structure:
            title = f"Brillouin Zone: {self.crystal_structure['name']} ({self.crystal_structure['lattice_type']})"
        else:
            title = "Brillouin Zone"
        
        ax.set_title(title)
        
        # Set limits
        if is_3d:
            ax.set_xlim([-1.2, 1.2])
            ax.set_ylim([-1.2, 1.2])
            ax.set_zlim([-1.2, 1.2])
        else:
            ax.set_xlim([-1.2, 1.2])
            ax.set_ylim([-1.2, 1.2])
        
        return fig
    
    def visualize_path(self, path, bands=None):
        """Visualize a path through the Brillouin zone.
        
        Args:
            path (list): List of special point names or coordinates
            bands (numpy.ndarray, optional): Band structure data along path
            
        Returns:
            tuple: (path_fig, bands_fig) - Two matplotlib figures
        """
        if not self.is_initialized() or not path:
            return None, None
        
        # Create figure for path visualization
        path_fig = Figure(figsize=(6, 5))
        path_ax = path_fig.add_subplot(111, projection='3d')
        
        # Get coordinates for path
        path_coords = []
        for point in path:
            if isinstance(point, str) and point in self.special_points:
                path_coords.append(self.special_points[point])
            elif isinstance(point, (list, tuple)) and len(point) == 3:
                path_coords.append(point)
        
        if not path_coords:
            return None, None
        
        # Convert to numpy array
        path_coords = np.array(path_coords)
        
        # Draw Brillouin zone
        vertices = self.brillouin_zone['vertices']
        edges = self.brillouin_zone['edges']
        
        for edge in edges:
            path_ax.plot3D(
                [vertices[edge[0], 0], vertices[edge[1], 0]],
                [vertices[edge[0], 1], vertices[edge[1], 1]],
                [vertices[edge[0], 2], vertices[edge[1], 2]],
                'b-', alpha=0.3
            )
        
        # Draw path
        path_ax.plot3D(
            path_coords[:, 0],
            path_coords[:, 1],
            path_coords[:, 2],
            'r-', linewidth=2
        )
        
        # Draw points along path
        path_ax.scatter(
            path_coords[:, 0],
            path_coords[:, 1],
            path_coords[:, 2],
            color='red', s=50
        )
        
        # Label points
        for i, point in enumerate(path):
            if isinstance(point, str):
                label = point
            else:
                label = f"P{i}"
            
            path_ax.text(
                path_coords[i, 0],
                path_coords[i, 1],
                path_coords[i, 2],
                label,
                fontsize=10
            )
        
        # Set labels and title
        path_ax.set_xlabel('kx')
        path_ax.set_ylabel('ky')
        path_ax.set_zlabel('kz')
        path_ax.set_title('Path through Brillouin Zone')
        
        # Set equal aspect ratio
        path_ax.set_box_aspect([1, 1, 1])
        
        # Create figure for band structure if provided
        bands_fig = None
        if bands is not None:
            bands_fig = Figure(figsize=(6, 5))
            bands_ax = bands_fig.add_subplot(111)
            
            # Generate x coordinates for the path
            x = np.linspace(0, 1, bands.shape[0])
            
            # Plot bands
            for i in range(bands.shape[1]):
                bands_ax.plot(x, bands[:, i], 'k-')
            
            # Add vertical lines and labels at special points
            x_ticks = [0]
            x_labels = [path[0] if isinstance(path[0], str) else 'P0']
            
            for i in range(1, len(path)):
                x_ticks.append(i / (len(path) - 1))
                x_labels.append(path[i] if isinstance(path[i], str) else f'P{i}')
                bands_ax.axvline(x=i / (len(path) - 1), color='k', linestyle='--', alpha=0.5)
            
            # Set labels, title and ticks
            bands_ax.set_xticks(x_ticks)
            bands_ax.set_xticklabels(x_labels)
            bands_ax.set_ylabel('Energy (eV)')
            bands_ax.set_title('Band Structure')
            bands_ax.grid(True, alpha=0.3)
        
        return path_fig, bands_fig 