#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=no-name-in-module, import-error
from PyQt5.QtWidgets import (
    QWidget,
    QGridLayout,
    QFormLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QTabWidget,
    QDoubleSpinBox,
    QGroupBox,
    QRadioButton,
    QButtonGroup,
    QFileDialog,
    QMessageBox,
    QComboBox,
    QListWidget,
    QTableWidget,
    QTableWidgetItem,
    QHeaderView,
    QSpinBox,
    QSplitter,
)
from PyQt5.QtCore import Qt, pyqtSlot
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

# Add parent directory to path to allow imports from sibling packages
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from packages.gui.tabs.tab_interface import TabInterface
from packages.trajectory_planner.interface import TrajectoryPlanner


class MatplotlibCanvas(FigureCanvas):
    """Matplotlib canvas for embedding in Qt applications."""
    
    def __init__(self, width=6, height=5, dpi=100, projection='3d'):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        if projection:
            self.axes = self.fig.add_subplot(111, projection=projection)
        else:
            self.axes = self.fig.add_subplot(111)
        super(MatplotlibCanvas, self).__init__(self.fig)
        self.fig.tight_layout()


class TrajectoryPlannerTab(TabInterface):
    """Tab for planning RIXS measurement trajectories."""

    def __init__(self, main_window=None):
        # Create backend instance
        self.planner = TrajectoryPlanner()
        self.planner.initialize()

        # Initialize UI
        super().__init__(main_window)

        # Set window title
        self.setWindowTitle("Trajectory Planner")

    def init_ui(self):
        """Initialize UI components."""
        # Create main splitter layout
        self.splitter = QSplitter(Qt.Horizontal)
        self.layout.addWidget(self.splitter, 0, 0)  # Add to grid at (0,0)

        # Left panel - controls
        left_widget = QWidget()
        left_layout = QGridLayout(left_widget)  # Changed to QGridLayout
        self.splitter.addWidget(left_widget)

        # Right panel - visualization
        right_widget = QWidget()
        right_layout = QGridLayout(right_widget)  # Changed to QGridLayout
        self.splitter.addWidget(right_widget)

        # Set default size ratio
        self.splitter.setSizes([400, 600])

        # Create left panel components
        self.create_point_editor(left_layout)
        self.create_trajectory_controls(left_layout)

        # Create visualization area
        self.create_visualization_area(right_layout)

    def create_point_editor(self, parent_layout):
        """Create the editor for trajectory points."""
        points_group = QGroupBox("Trajectory Points")
        points_layout = QGridLayout(points_group)  # Changed to QGridLayout

        # Table for points
        self.points_table = QTableWidget()
        self.points_table.setColumnCount(5)
        self.points_table.setHorizontalHeaderLabels(["H", "K", "L", "Time (s)", "Temp (K)"])
        self.points_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        points_layout.addWidget(self.points_table, 0, 0)  # Add to grid at (0,0)

        # Buttons for point manipulation
        buttons_layout = QGridLayout()  # Changed to QGridLayout

        add_button = QPushButton("Add Point")
        add_button.clicked.connect(self.add_point)
        buttons_layout.addWidget(add_button, 0, 0)  # Add to grid at (0,0)

        remove_button = QPushButton("Remove Point")
        remove_button.clicked.connect(self.remove_point)
        buttons_layout.addWidget(remove_button, 0, 1)  # Add to grid at (0,1)

        clear_button = QPushButton("Clear All")
        clear_button.clicked.connect(self.clear_points)
        buttons_layout.addWidget(clear_button, 0, 2)  # Add to grid at (0,2)

        points_layout.addLayout(buttons_layout, 1, 0)  # Add to grid at (1,0)

        # Add point form
        form_group = QGroupBox("Add Point")
        form_layout = QFormLayout(form_group)

        self.h_input = QDoubleSpinBox()
        self.h_input.setRange(-10.0, 10.0)
        self.h_input.setDecimals(3)
        self.h_input.setValue(0.0)
        form_layout.addRow("H:", self.h_input)

        self.k_input = QDoubleSpinBox()
        self.k_input.setRange(-10.0, 10.0)
        self.k_input.setDecimals(3)
        self.k_input.setValue(0.0)
        form_layout.addRow("K:", self.k_input)

        self.l_input = QDoubleSpinBox()
        self.l_input.setRange(-10.0, 10.0)
        self.l_input.setDecimals(3)
        self.l_input.setValue(0.0)
        form_layout.addRow("L:", self.l_input)

        self.time_input = QDoubleSpinBox()
        self.time_input.setRange(1.0, 3600.0)
        self.time_input.setValue(60.0)
        self.time_input.setSuffix(" s")
        form_layout.addRow("Integration time:", self.time_input)

        self.temp_input = QDoubleSpinBox()
        self.temp_input.setRange(0.0, 1000.0)
        self.temp_input.setValue(300.0)
        self.temp_input.setSuffix(" K")
        form_layout.addRow("Temperature:", self.temp_input)

        points_layout.addWidget(form_group, 2, 0)  # Add to grid at (2,0)

        parent_layout.addWidget(points_group, 0, 0)  # Add to grid at (0,0)

    def create_trajectory_controls(self, parent_layout):
        """Create controls for trajectory calculation and export."""
        trajectory_group = QGroupBox("Trajectory Settings")
        trajectory_layout = QGridLayout(trajectory_group)  # Changed to QGridLayout

        # Temperature profile form
        temp_form = QFormLayout()

        self.start_temp_input = QDoubleSpinBox()
        self.start_temp_input.setRange(0.0, 1000.0)
        self.start_temp_input.setValue(300.0)
        self.start_temp_input.setSuffix(" K")
        temp_form.addRow("Start temperature:", self.start_temp_input)

        self.end_temp_input = QDoubleSpinBox()
        self.end_temp_input.setRange(0.0, 1000.0)
        self.end_temp_input.setValue(300.0)
        self.end_temp_input.setSuffix(" K")
        temp_form.addRow("End temperature:", self.end_temp_input)

        self.temp_steps_input = QSpinBox()
        self.temp_steps_input.setRange(2, 100)
        self.temp_steps_input.setValue(10)
        temp_form.addRow("Temperature steps:", self.temp_steps_input)

        trajectory_layout.addLayout(temp_form, 0, 0)  # Add to grid at (0,0)

        # Trajectory mode
        mode_form = QFormLayout()
        self.trajectory_mode_combo = QComboBox()
        self.trajectory_mode_combo.addItems(["Linear", "Zigzag", "Spiral"])
        mode_form.addRow("Trajectory mode:", self.trajectory_mode_combo)

        trajectory_layout.addLayout(mode_form, 1, 0)  # Add to grid at (1,0)

        # Calculate button
        calculate_button = QPushButton("Calculate Trajectory")
        calculate_button.clicked.connect(self.calculate_trajectory)
        trajectory_layout.addWidget(calculate_button, 2, 0)  # Add to grid at (2,0)

        # Export form
        export_form = QFormLayout()
        self.export_format_combo = QComboBox()
        self.export_format_combo.addItems(["JSON", "SPEC", "Bluesky"])
        export_form.addRow("Export format:", self.export_format_combo)

        export_button = QPushButton("Export Trajectory")
        export_button.clicked.connect(self.export_trajectory)

        trajectory_layout.addLayout(export_form, 3, 0)  # Add to grid at (3,0)
        trajectory_layout.addWidget(export_button, 4, 0)  # Add to grid at (4,0)

        parent_layout.addWidget(trajectory_group, 1, 0)  # Add to grid at (1,0)

    def create_visualization_area(self, parent_layout):
        """Create area for visualizing the trajectory."""
        vis_group = QGroupBox("Visualization")
        vis_layout = QGridLayout(vis_group)  # Changed to QGridLayout

        # Add matplotlib canvas
        self.canvas = MatplotlibCanvas(width=8, height=6)
        vis_layout.addWidget(self.canvas, 0, 0)  # Add to grid at (0,0)

        parent_layout.addWidget(vis_group, 0, 0)  # Add to grid at (0,0)

    @pyqtSlot()
    def add_point(self):
        """Add a point to the trajectory."""
        h = self.h_input.value()
        k = self.k_input.value()
        l = self.l_input.value()
        time = self.time_input.value()
        temp = self.temp_input.value()

        try:
            index = self.planner.add_point(h, k, l, time, temp)
            self.update_points_table()
            self.update_visualizations()

            # Give feedback
            self.statusBar().showMessage(f"Added point ({h}, {k}, {l})", 2000)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error adding point: {str(e)}")

    @pyqtSlot()
    def remove_point(self):
        """Remove selected point from trajectory."""
        selected = self.points_table.selectedIndexes()
        if not selected:
            return

        # Get row index (take first selected item)
        row = selected[0].row()

        try:
            if self.planner.remove_point(row):
                self.update_points_table()
                self.update_visualizations()
                self.statusBar().showMessage("Point removed", 2000)
            else:
                QMessageBox.warning(self, "Warning", "Could not remove point")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error removing point: {str(e)}")

    @pyqtSlot()
    def clear_points(self):
        """Clear all trajectory points."""
        reply = QMessageBox.question(
            self, "Confirm Clear", "Are you sure you want to clear all points?",
            QMessageBox.Yes | QMessageBox.No, QMessageBox.No
        )

        if reply == QMessageBox.Yes:
            self.planner.initialize()  # Reset the planner
            self.update_points_table()
            self.update_visualizations()
            self.statusBar().showMessage("All points cleared", 2000)

    def update_points_table(self):
        """Update the table of trajectory points."""
        points = self.planner.get_points()

        # Update table
        self.points_table.setRowCount(len(points))

        for i, point in enumerate(points):
            # H, K, L
            self.points_table.setItem(i, 0, QTableWidgetItem(f"{point['h']:.3f}"))
            self.points_table.setItem(i, 1, QTableWidgetItem(f"{point['k']:.3f}"))
            self.points_table.setItem(i, 2, QTableWidgetItem(f"{point['l']:.3f}"))

            # Integration time
            time_item = QTableWidgetItem(f"{point.get('integration_time', 60.0):.1f}")
            self.points_table.setItem(i, 3, time_item)

            # Temperature (may be None)
            temp = point.get('temperature')
            temp_item = QTableWidgetItem(f"{temp:.1f}" if temp is not None else "")
            self.points_table.setItem(i, 4, temp_item)

    @pyqtSlot()
    def calculate_trajectory(self):
        """Calculate the trajectory based on current settings."""
        if not self.planner.get_points():
            QMessageBox.warning(self, "Warning", "No trajectory points defined")
            return

        try:
            # Get temperature settings
            start_temp = self.start_temp_input.value()
            end_temp = self.end_temp_input.value()
            steps = self.temp_steps_input.value()

            # Get trajectory mode
            mode_text = self.trajectory_mode_combo.currentText().lower()

            # Calculate trajectory
            result = self.planner.calculate_trajectory(
                start_temp=start_temp,
                end_temp=end_temp,
                steps=steps,
                mode=mode_text
            )

            if result:
                # Update visualizations
                self.update_visualizations()

                # Show summary
                msg = (f"Trajectory calculated with {len(result['points'])} points\n"
                      f"Total distance: {result['total_distance']:.2f}\n"
                      f"Estimated time: {result['total_time']:.1f} seconds")
                QMessageBox.information(self, "Trajectory Calculated", msg)
            else:
                QMessageBox.warning(self, "Warning", "Failed to calculate trajectory")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error calculating trajectory: {str(e)}")

    @pyqtSlot()
    def export_trajectory(self):
        """Export the trajectory in the selected format."""
        if not self.planner.get_points():
            QMessageBox.warning(self, "Warning", "No trajectory points defined")
            return

        try:
            # Get export format
            format_text = self.export_format_combo.currentText().lower()

            # Get save path
            file_types = {
                "json": "JSON Files (*.json)",
                "spec": "SPEC Files (*.mac)",
                "bluesky": "Python Files (*.py)"
            }

            file_path, _ = QFileDialog.getSaveFileName(
                self, "Export Trajectory", "", file_types.get(format_text, "All Files (*)")
            )

            if not file_path:
                return

            # Export trajectory
            exported = self.planner.export_trajectory(format_type=format_text)

            if exported:
                # For this template, just show the result
                # In a real implementation, we would write to the file
                if format_text == "spec":
                    # Write text to file
                    with open(file_path, 'w') as f:
                        f.write(exported)
                    QMessageBox.information(self, "Export Successful", 
                                          f"Trajectory exported to {file_path}")
                else:
                    # For JSON or other formats, we would serialize the data
                    QMessageBox.information(self, "Export Successful", 
                                          "Trajectory export format not fully implemented")
            else:
                QMessageBox.warning(self, "Warning", "Failed to export trajectory")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error exporting trajectory: {str(e)}")

    def update_visualizations(self):
        """Update all visualizations."""
        self.update_3d_visualization()
        self.update_2d_visualization()
        self.update_temperature_visualization()

    def update_3d_visualization(self):
        """Update the 3D visualization."""
        if not self.planner.get_points():
            # Clear the canvas
            self.canvas.axes.clear()
            self.canvas.draw()
            return

        # Generate figure with 3D view
        fig = self.planner.visualize_trajectory(view_direction=None)

        if fig:
            # Clear current axes
            self.canvas.axes.clear()

            # Copy content from generated figure to canvas figure
            for ax in fig.get_axes():
                # Copy collections (points)
                for collection in ax.collections:
                    self.canvas.axes.add_collection3d(collection.copy())

                # Copy lines
                for line in ax.lines:
                    self.canvas.axes.plot3D(
                        line.get_xdata(),
                        line.get_ydata(),
                        line.get_zdata(),
                        color=line.get_color(),
                        linestyle=line.get_linestyle(),
                        linewidth=line.get_linewidth(),
                    )

                # Copy texts
                for text in ax.texts:
                    self.canvas.axes.text(
                        text.get_position()[0],
                        text.get_position()[1],
                        text.get_position()[2],
                        text.get_text(),
                        fontsize=text.get_fontsize(),
                    )

                # Copy axis settings
                self.canvas.axes.set_xlim(ax.get_xlim())
                self.canvas.axes.set_ylim(ax.get_ylim())
                self.canvas.axes.set_zlim(ax.get_zlim())
                self.canvas.axes.set_xlabel(ax.get_xlabel())
                self.canvas.axes.set_ylabel(ax.get_ylabel())
                self.canvas.axes.set_zlabel(ax.get_zlabel())
                self.canvas.axes.set_title(ax.get_title())

            self.canvas.draw()

    @pyqtSlot(int)
    def update_2d_visualization(self, index=None):
        """Update the 2D visualization."""
        if not self.planner.get_points():
            # Clear the canvas
            self.canvas.axes.clear()
            self.canvas.draw()
            return

        # Get view direction based on selected index
        view_index = self.view_direction_combo.currentIndex()
        if view_index == 0:
            view_direction = 'z'  # HK plane
        elif view_index == 1:
            view_direction = 'y'  # HL plane
        else:
            view_direction = 'x'  # KL plane

        # Generate figure with 2D view
        fig = self.planner.visualize_trajectory(view_direction=view_direction)

        if fig:
            # Clear current axes
            self.canvas.axes.clear()

            # Copy content from generated figure to canvas figure
            for ax in fig.get_axes():
                # Copy collections (points)
                for collection in ax.collections:
                    self.canvas.axes.add_collection(collection.copy())

                # Copy lines
                for line in ax.lines:
                    self.canvas.axes.plot(
                        line.get_xdata(),
                        line.get_ydata(),
                        color=line.get_color(),
                        linestyle=line.get_linestyle(),
                        linewidth=line.get_linewidth(),
                    )

                # Copy texts
                for text in ax.texts:
                    self.canvas.axes.text(
                        text.get_position()[0],
                        text.get_position()[1],
                        text.get_text(),
                        fontsize=text.get_fontsize(),
                    )

                # Copy axis settings
                self.canvas.axes.set_xlim(ax.get_xlim())
                self.canvas.axes.set_ylim(ax.get_ylim())
                self.canvas.axes.set_xlabel(ax.get_xlabel())
                self.canvas.axes.set_ylabel(ax.get_ylabel())
                self.canvas.axes.set_title(ax.get_title())

                # Set aspect ratio
                if hasattr(ax, 'get_aspect') and ax.get_aspect() == 'equal':
                    self.canvas.axes.set_aspect("equal")

            self.canvas.draw()

    def update_temperature_visualization(self):
        """Update the temperature profile visualization."""
        # Clear current axes
        self.canvas.axes.clear()

        # Check if we have a calculated trajectory with temperature profile
        if not hasattr(self.planner, 'current_trajectory') or not self.planner.current_trajectory:
            self.canvas.draw()
            return

        # Generate temperature profile figure
        fig = self.planner.visualize_temperature_profile()

        if fig:
            # Copy content from generated figure to canvas figure
            for ax in fig.get_axes():
                # Copy collections (points)
                for collection in ax.collections:
                    self.canvas.axes.add_collection(collection.copy())

                # Copy lines
                for line in ax.lines:
                    self.canvas.axes.plot(
                        line.get_xdata(),
                        line.get_ydata(),
                        color=line.get_color(),
                        linestyle=line.get_linestyle(),
                        linewidth=line.get_linewidth(),
                    )

                # Copy axis settings
                self.canvas.axes.set_xlim(ax.get_xlim())
                self.canvas.axes.set_ylim(ax.get_ylim())
                self.canvas.axes.set_xlabel(ax.get_xlabel())
                self.canvas.axes.set_ylabel(ax.get_ylabel())
                self.canvas.axes.set_title(ax.get_title())

                # Copy grid
                self.canvas.axes.grid(True, alpha=0.3)

        self.canvas.draw()

    def statusBar(self):
        """Get the main window's status bar."""
        try:
            return self.parent().parent().parent().statusBar()
        except:
            # Fallback if we can't get the main window's status bar
            class DummyStatusBar:
                def showMessage(self, msg, timeout=0):
                    pass
            return DummyStatusBar()

    def get_module_instance(self):
        """Get the backend module instance."""
        return self.planner 
