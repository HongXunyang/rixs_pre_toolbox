#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from matplotlib.figure import Figure


class TrajectoryPlanner:
    """Interface for the RIXS trajectory planner.
    
    This class handles all the calculations required for planning measurement
    trajectories. It's a pure Python implementation without PyQt dependencies.
    """
    
    def __init__(self):
        """Initialize the planner."""
        self._initialized = False
        self.trajectory_points = []
        self.current_trajectory = None
        self.metadata = {}
    
    def initialize(self, points=None, metadata=None):
        """Initialize the planner with initial trajectory points.
        
        Args:
            points (list, optional): Initial trajectory points
            metadata (dict, optional): Additional metadata for the trajectory
            
        Returns:
            bool: True if initialization was successful
        """
        try:
            self.trajectory_points = points or []
            self.metadata = metadata or {}
            self._initialized = True
            return True
        except Exception as e:
            print(f"Error initializing trajectory planner: {str(e)}")
            return False
    
    def is_initialized(self):
        """Check if the planner is initialized.
        
        Returns:
            bool: True if the planner is initialized
        """
        return self._initialized
    
    def add_point(self, h, k, l, integration_time=None, temperature=None, custom_params=None):
        """Add a point to the trajectory.
        
        Args:
            h, k, l (float): HKL coordinates of the point
            integration_time (float, optional): Integration time in seconds
            temperature (float, optional): Temperature in Kelvin
            custom_params (dict, optional): Custom parameters for this point
            
        Returns:
            int: Index of the added point
        """
        if not self.is_initialized():
            raise ValueError("Planner not initialized")
        
        point = {
            'h': h,
            'k': k,
            'l': l,
            'integration_time': integration_time or 60.0,  # Default 60 seconds
            'temperature': temperature,
            'custom_params': custom_params or {}
        }
        
        self.trajectory_points.append(point)
        return len(self.trajectory_points) - 1
    
    def remove_point(self, index):
        """Remove a point from the trajectory.
        
        Args:
            index (int): Index of the point to remove
            
        Returns:
            bool: True if the point was removed successfully
        """
        if not self.is_initialized() or index < 0 or index >= len(self.trajectory_points):
            return False
        
        self.trajectory_points.pop(index)
        return True
    
    def update_point(self, index, **kwargs):
        """Update properties of a trajectory point.
        
        Args:
            index (int): Index of the point to update
            **kwargs: Properties to update
            
        Returns:
            bool: True if the point was updated successfully
        """
        if not self.is_initialized() or index < 0 or index >= len(self.trajectory_points):
            return False
        
        point = self.trajectory_points[index]
        for key, value in kwargs.items():
            if key in point:
                point[key] = value
        
        return True
    
    def get_points(self):
        """Get all trajectory points.
        
        Returns:
            list: List of all trajectory points
        """
        return self.trajectory_points
    
    def calculate_trajectory(self, start_temp=None, end_temp=None, steps=None, mode="linear"):
        """Calculate a trajectory through the defined points.
        
        Args:
            start_temp (float, optional): Starting temperature
            end_temp (float, optional): Ending temperature
            steps (int, optional): Number of temperature steps
            mode (str): Trajectory mode ('linear', 'zigzag', 'spiral')
            
        Returns:
            dict: Calculated trajectory information
        """
        if not self.is_initialized() or not self.trajectory_points:
            return None
        
        # Create a simple trajectory
        # In a real implementation, this would calculate optimal paths, 
        # motor movements, etc.
        
        # Simple linear interpolation between points
        points = np.array([[p['h'], p['k'], p['l']] for p in self.trajectory_points])
        
        # Calculate total distance
        distances = np.zeros(len(points))
        for i in range(1, len(points)):
            distances[i] = distances[i-1] + np.linalg.norm(points[i] - points[i-1])
        
        total_distance = distances[-1]
        
        # Generate temperature steps if provided
        temp_profile = None
        if start_temp is not None and end_temp is not None and steps is not None:
            if mode == "linear":
                temp_profile = np.linspace(start_temp, end_temp, steps)
            else:
                # Other temperature profiles could be implemented
                temp_profile = np.linspace(start_temp, end_temp, steps)
        
        # Calculate estimated time
        total_time = sum(p.get('integration_time', 60.0) for p in self.trajectory_points)
        
        # Store the calculated trajectory
        self.current_trajectory = {
            'points': self.trajectory_points,
            'distances': distances.tolist(),
            'total_distance': total_distance,
            'total_time': total_time,
            'temperature_profile': temp_profile.tolist() if temp_profile is not None else None,
            'mode': mode
        }
        
        return self.current_trajectory
    
    def export_trajectory(self, format_type="json"):
        """Export the trajectory in the specified format.
        
        Args:
            format_type (str): Export format ('json', 'spec', 'bluesky')
            
        Returns:
            str: Trajectory in the specified format
        """
        if not self.is_initialized() or not self.current_trajectory:
            return None
        
        # This is a placeholder implementation
        # In a real implementation, this would generate the appropriate
        # format for the specified export type
        
        if format_type == "json":
            # Return a simplified JSON representation
            return {
                'points': self.trajectory_points,
                'metadata': self.metadata
            }
        elif format_type == "spec":
            # Return a SPEC macro format (simplified)
            lines = ["# RIXS trajectory generated by TrajectoryPlanner"]
            lines.append(f"# Total points: {len(self.trajectory_points)}")
            lines.append("# H K L Time Temp")
            
            for p in self.trajectory_points:
                lines.append(f"mv_hkl({p['h']}, {p['k']}, {p['l']})")
                lines.append(f"count {p.get('integration_time', 60.0)}")
            
            return "\n".join(lines)
        else:
            # Default to JSON
            return {
                'points': self.trajectory_points,
                'metadata': self.metadata
            }
    
    def visualize_trajectory(self, view_direction=None):
        """Visualize the trajectory.
        
        Args:
            view_direction (str): View direction ('x', 'y', 'z' or None for 3D)
            
        Returns:
            Figure: Matplotlib figure with the visualization
        """
        if not self.is_initialized() or not self.trajectory_points:
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
        
        # Extract point coordinates
        points = np.array([[p['h'], p['k'], p['l']] for p in self.trajectory_points])
        
        if is_3d:
            # Draw trajectory in 3D
            ax.plot3D(points[:, 0], points[:, 1], points[:, 2], 'r-', linewidth=2)
            ax.scatter(points[:, 0], points[:, 1], points[:, 2], color='blue', s=50)
            
            # Label points
            for i, point in enumerate(points):
                ax.text(point[0], point[1], point[2], str(i), fontsize=10)
            
            # Set labels
            ax.set_xlabel('H')
            ax.set_ylabel('K')
            ax.set_zlabel('L')
            
            # Set equal aspect ratio
            ax.set_box_aspect([1, 1, 1])
        else:
            # Create 2D projection
            if view_direction == 'x':
                idx1, idx2 = 1, 2  # K, L
                xlabel, ylabel = 'K', 'L'
            elif view_direction == 'y':
                idx1, idx2 = 0, 2  # H, L
                xlabel, ylabel = 'H', 'L'
            else:  # 'z' or default
                idx1, idx2 = 0, 1  # H, K
                xlabel, ylabel = 'H', 'K'
            
            # Draw projection of trajectory
            ax.plot(points[:, idx1], points[:, idx2], 'r-', linewidth=2)
            ax.scatter(points[:, idx1], points[:, idx2], color='blue', s=50)
            
            # Label points
            for i, point in enumerate(points):
                ax.text(point[idx1], point[idx2], str(i), fontsize=10)
            
            # Set labels
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            
            # Set equal aspect ratio
            ax.set_aspect('equal')
        
        # Set title
        ax.set_title('Measurement Trajectory')
        
        return fig
    
    def visualize_temperature_profile(self):
        """Visualize the temperature profile if available.
        
        Returns:
            Figure: Matplotlib figure with the temperature profile
        """
        if (not self.is_initialized() or 
            not self.current_trajectory or 
            self.current_trajectory.get('temperature_profile') is None):
            return None
        
        # Create figure
        fig = Figure(figsize=(6, 4))
        ax = fig.add_subplot(111)
        
        # Get temperature profile
        temps = self.current_trajectory['temperature_profile']
        x = np.linspace(0, 1, len(temps))
        
        # Plot temperature profile
        ax.plot(x, temps, 'r-', linewidth=2)
        ax.scatter(x, temps, color='blue', s=30)
        
        # Set labels and title
        ax.set_xlabel('Normalized Trajectory Progress')
        ax.set_ylabel('Temperature (K)')
        ax.set_title('Temperature Profile')
        
        # Add grid
        ax.grid(True, alpha=0.3)
        
        return fig 