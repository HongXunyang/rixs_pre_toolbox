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
    QTextEdit,
    QSplitter,
    QSpinBox,
    QCheckBox,
)
from PyQt5.QtCore import Qt, pyqtSlot
import sys
import os
import numpy as np
import random
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

# Add parent directory to path to allow imports from sibling packages
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from packages.gui.tabs.tab_interface import TabInterface
from packages.placeholder_module.interface import PlaceholderModule


class MatplotlibCanvas(FigureCanvas):
    """Matplotlib canvas for embedding in Qt applications."""
    
    def __init__(self, width=6, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        super(MatplotlibCanvas, self).__init__(self.fig)
        self.fig.tight_layout()


class PlaceholderTab(TabInterface):
    """Tab for demonstrating how to add new tabs to the application."""

    def __init__(self, main_window=None):
        # Create backend instance
        self.module = PlaceholderModule()
        self.module.initialize()

        # Initialize UI
        super().__init__(main_window)

        # Set window title
        self.setWindowTitle("Placeholder Module")

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

        # Create panels
        self.create_input_panel(left_layout)
        self.create_processing_panel(left_layout)
        self.create_visualization_panel(right_layout)

    def create_input_panel(self, parent_layout):
        """Create the input panel for the tab."""
        input_group = QGroupBox("Input Data")
        input_layout = QGridLayout(input_group)  # Changed to QGridLayout

        # Data type selection
        data_type_form = QFormLayout()
        self.data_type_combo = QComboBox()
        self.data_type_combo.addItems(["Random Array", "Sine Wave", "Custom Input"])
        self.data_type_combo.currentIndexChanged.connect(self.update_data_form)
        data_type_form.addRow("Data Type:", self.data_type_combo)
        input_layout.addLayout(data_type_form, 0, 0)  # Add to grid at (0,0)

        # Data parameters group
        self.data_params_group = QGroupBox("Data Parameters")
        self.data_params_layout = QFormLayout(self.data_params_group)

        # Random array parameters (initial view)
        self.size_input = QSpinBox()
        self.size_input.setRange(10, 1000)
        self.size_input.setValue(100)
        self.data_params_layout.addRow("Size:", self.size_input)

        self.min_input = QDoubleSpinBox()
        self.min_input.setRange(-1000, 1000)
        self.min_input.setValue(0)
        self.data_params_layout.addRow("Min Value:", self.min_input)

        self.max_input = QDoubleSpinBox()
        self.max_input.setRange(-1000, 1000)
        self.max_input.setValue(100)
        self.data_params_layout.addRow("Max Value:", self.max_input)

        input_layout.addWidget(self.data_params_group, 1, 0)  # Add to grid at (1,0)

        # Custom input for text data (initially hidden)
        self.custom_input_group = QGroupBox("Custom Input")
        custom_input_layout = QGridLayout(
            self.custom_input_group
        )  # Changed to QGridLayout

        self.custom_input_text = QTextEdit()
        self.custom_input_text.setPlaceholderText("Enter custom data here (comma-separated values)")
        custom_input_layout.addWidget(
            self.custom_input_text, 0, 0
        )  # Add to grid at (0,0)

        input_layout.addWidget(self.custom_input_group, 2, 0)  # Add to grid at (2,0)
        self.custom_input_group.setVisible(False)

        # Generate button
        generate_button = QPushButton("Generate Data")
        generate_button.clicked.connect(self.generate_data)
        input_layout.addWidget(generate_button, 3, 0)  # Add to grid at (3,0)

        # Add to parent layout
        parent_layout.addWidget(input_group, 0, 0)  # Add to grid at (0,0)

        # Set initial form state
        self.update_data_form(0)  # 0 = Random Array

    def create_processing_panel(self, parent_layout):
        """Create the processing panel for the tab."""
        processing_group = QGroupBox("Processing")
        processing_layout = QGridLayout(processing_group)  # Changed to QGridLayout

        # Processing type selection
        process_type_form = QFormLayout()
        self.process_type_combo = QComboBox()
        self.process_type_combo.addItems(["None", "Smooth", "Filter", "Transform"])
        process_type_form.addRow("Process Type:", self.process_type_combo)
        processing_layout.addLayout(process_type_form, 0, 0)  # Add to grid at (0,0)

        # Processing parameters
        params_form = QFormLayout()

        self.param1_input = QDoubleSpinBox()
        self.param1_input.setRange(0, 100)
        self.param1_input.setValue(1.0)
        params_form.addRow("Parameter 1:", self.param1_input)

        self.param2_input = QDoubleSpinBox()
        self.param2_input.setRange(0, 100)
        self.param2_input.setValue(1.0)
        params_form.addRow("Parameter 2:", self.param2_input)

        processing_layout.addLayout(params_form, 1, 0)  # Add to grid at (1,0)

        # Additional options
        self.normalize_check = QCheckBox("Normalize Output")
        processing_layout.addWidget(self.normalize_check, 2, 0)  # Add to grid at (2,0)

        # Process button
        process_button = QPushButton("Process Data")
        process_button.clicked.connect(self.process_data)
        processing_layout.addWidget(process_button, 3, 0)  # Add to grid at (3,0)

        # Add to parent layout
        parent_layout.addWidget(processing_group, 1, 0)  # Add to grid at (1,0)

    def create_visualization_panel(self, parent_layout):
        """Create the visualization panel for the tab."""
        visualization_group = QGroupBox("Visualization")
        visualization_layout = QGridLayout(
            visualization_group
        )  # Changed to QGridLayout

        # Visualization type selection
        vis_type_form = QFormLayout()
        self.vis_type_combo = QComboBox()
        self.vis_type_combo.addItems(["Line Plot", "Bar Chart", "Scatter Plot"])
        self.vis_type_combo.currentIndexChanged.connect(self.update_visualization)
        vis_type_form.addRow("Chart Type:", self.vis_type_combo)
        visualization_layout.addLayout(vis_type_form, 0, 0)  # Add to grid at (0,0)

        # Plot customization
        customize_form = QFormLayout()

        self.title_input = QLineEdit("Placeholder Data")
        customize_form.addRow("Title:", self.title_input)

        self.color_combo = QComboBox()
        self.color_combo.addItems(["Blue", "Red", "Green", "Orange", "Purple"])
        customize_form.addRow("Color:", self.color_combo)

        self.grid_check = QCheckBox("Show Grid")
        self.grid_check.setChecked(True)
        customize_form.addRow("", self.grid_check)

        visualization_layout.addLayout(customize_form, 1, 0)  # Add to grid at (1,0)

        # Matplotlib canvas for visualization
        self.canvas = MatplotlibCanvas()
        visualization_layout.addWidget(self.canvas, 2, 0)  # Add to grid at (2,0)

        # Export button
        export_button = QPushButton("Export Plot")
        export_button.clicked.connect(self.export_plot)
        visualization_layout.addWidget(export_button, 3, 0)  # Add to grid at (3,0)

        # Add to parent layout
        parent_layout.addWidget(visualization_group, 0, 0)  # Add to grid at (0,0)

    @pyqtSlot(int)
    def update_data_form(self, index):
        """Update the data form based on the selected data type."""
        # Hide both forms initially
        self.data_params_group.setVisible(False)
        self.custom_input_group.setVisible(False)

        # Show the appropriate form
        if index == 0:  # Random Array
            self.setup_random_array_form()
            self.data_params_group.setVisible(True)
        elif index == 1:  # Sine Wave
            self.setup_sine_wave_form()
            self.data_params_group.setVisible(True)
        elif index == 2:  # Custom Input
            self.custom_input_group.setVisible(True)

    def setup_random_array_form(self):
        """Set up the form for random array data."""
        # Clear existing layout
        self.clear_layout(self.data_params_layout)

        # Add random array parameters
        self.size_input = QSpinBox()
        self.size_input.setRange(10, 1000)
        self.size_input.setValue(100)
        self.data_params_layout.addRow("Size:", self.size_input)

        self.min_input = QDoubleSpinBox()
        self.min_input.setRange(-1000, 1000)
        self.min_input.setValue(0)
        self.data_params_layout.addRow("Min Value:", self.min_input)

        self.max_input = QDoubleSpinBox()
        self.max_input.setRange(-1000, 1000)
        self.max_input.setValue(100)
        self.data_params_layout.addRow("Max Value:", self.max_input)

    def setup_sine_wave_form(self):
        """Set up the form for sine wave data."""
        # Clear existing layout
        self.clear_layout(self.data_params_layout)

        # Add sine wave parameters
        self.size_input = QSpinBox()
        self.size_input.setRange(10, 1000)
        self.size_input.setValue(100)
        self.data_params_layout.addRow("Points:", self.size_input)

        self.amplitude_input = QDoubleSpinBox()
        self.amplitude_input.setRange(0.1, 100)
        self.amplitude_input.setValue(1.0)
        self.data_params_layout.addRow("Amplitude:", self.amplitude_input)

        self.frequency_input = QDoubleSpinBox()
        self.frequency_input.setRange(0.1, 50)
        self.frequency_input.setValue(1.0)
        self.frequency_input.setSingleStep(0.1)
        self.data_params_layout.addRow("Frequency:", self.frequency_input)

        self.noise_input = QDoubleSpinBox()
        self.noise_input.setRange(0, 1)
        self.noise_input.setValue(0.1)
        self.noise_input.setSingleStep(0.05)
        self.data_params_layout.addRow("Noise:", self.noise_input)

    def clear_layout(self, layout):
        """Clear all items from a layout."""
        if layout is None:
            return

        while layout.count():
            item = layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()
            else:
                self.clear_layout(item.layout())

    @pyqtSlot()
    def generate_data(self):
        """Generate data based on the selected type."""
        try:
            data_type = self.data_type_combo.currentIndex()

            if data_type == 0:  # Random Array
                size = self.size_input.value()
                min_val = self.min_input.value()
                max_val = self.max_input.value()

                data = np.random.uniform(min_val, max_val, size)
                self.module.set_data(data)

                self.statusBar().showMessage(f"Generated random array with {size} points", 2000)

            elif data_type == 1:  # Sine Wave
                points = self.size_input.value()
                amplitude = self.amplitude_input.value()
                frequency = self.frequency_input.value()
                noise = self.noise_input.value()

                x = np.linspace(0, 2 * np.pi, points)
                data = amplitude * np.sin(frequency * x)

                # Add noise if specified
                if noise > 0:
                    data += np.random.normal(0, noise * amplitude, points)

                self.module.set_data(data)

                self.statusBar().showMessage(f"Generated sine wave with {points} points", 2000)

            elif data_type == 2:  # Custom Input
                text = self.custom_input_text.toPlainText()
                if not text:
                    QMessageBox.warning(self, "Warning", "Please enter some data")
                    return

                try:
                    # Try to parse as comma-separated values
                    values = [float(x.strip()) for x in text.split(',') if x.strip()]
                    if not values:
                        raise ValueError("No valid values found")

                    data = np.array(values)
                    self.module.set_data(data)

                    self.statusBar().showMessage(f"Parsed {len(data)} custom values", 2000)
                except ValueError as e:
                    QMessageBox.warning(self, "Error", f"Could not parse input: {str(e)}")
                    return

            # Update visualization
            self.update_visualization()

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error generating data: {str(e)}")

    @pyqtSlot()
    def process_data(self):
        """Process the data with the selected method."""
        try:
            # Check if we have data
            if self.module.get_data() is None:
                QMessageBox.warning(self, "Warning", "No data to process")
                return

            process_type = self.process_type_combo.currentText()
            param1 = self.param1_input.value()
            param2 = self.param2_input.value()
            normalize = self.normalize_check.isChecked()

            # Process data based on selected type
            if process_type == "None":
                # No processing, just return
                return

            # Process the data
            result = self.module.process(
                process_type=process_type,
                param1=param1,
                param2=param2,
                normalize=normalize
            )

            if result:
                # In a real implementation, this would update the data in the module
                # For this template, we'll just show a message
                QMessageBox.information(
                    self, 
                    "Processing Complete", 
                    f"Data processed with {process_type}\n"
                    f"Parameters: {param1}, {param2}\n"
                    f"Normalize: {normalize}"
                )

                # Update visualization
                self.update_visualization()

                self.statusBar().showMessage(f"Data processed with {process_type}", 2000)
            else:
                QMessageBox.warning(self, "Warning", "Failed to process data")

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error processing data: {str(e)}")

    @pyqtSlot(int)
    def update_visualization(self, index=None):
        """Update the visualization with current data and settings."""
        try:
            # Check if we have data
            data = self.module.get_data()
            if data is None:
                # Generate placeholder visualization with no data
                fig = self.module.visualize()

                # Clear current axes
                self.canvas.axes.clear()

                # Copy text elements from the figure
                for ax in fig.get_axes():
                    for text in ax.texts:
                        self.canvas.axes.text(
                            text.get_position()[0], 
                            text.get_position()[1], 
                            text.get_text(),
                            fontsize=text.get_fontsize(),
                            ha='center', va='center'
                        )

                    self.canvas.axes.set_xlim(ax.get_xlim())
                    self.canvas.axes.set_ylim(ax.get_ylim())
                    self.canvas.axes.set_title(ax.get_title())
                    self.canvas.axes.axis('off')

                self.canvas.draw()
                return

            # Clear existing plot
            self.canvas.axes.clear()

            # Get visualization settings
            vis_type = self.vis_type_combo.currentText()
            title = self.title_input.text()
            color = self.color_combo.currentText().lower()
            show_grid = self.grid_check.isChecked()

            # Generate x-axis values
            x = np.arange(len(data))

            # Create plot based on type
            if vis_type == "Line Plot":
                self.canvas.axes.plot(x, data, f"{color}-", linewidth=2)
            elif vis_type == "Bar Chart":
                self.canvas.axes.bar(x, data, color=color, alpha=0.7)
            elif vis_type == "Scatter Plot":
                self.canvas.axes.scatter(x, data, color=color, alpha=0.7)

            # Set labels and title
            self.canvas.axes.set_xlabel("Index")
            self.canvas.axes.set_ylabel("Value")
            self.canvas.axes.set_title(title)

            # Set grid
            self.canvas.axes.grid(show_grid, alpha=0.3)

            # Update the canvas
            self.canvas.fig.tight_layout()
            self.canvas.draw()

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error updating visualization: {str(e)}")

    @pyqtSlot()
    def export_plot(self):
        """Export the current plot to a file."""
        try:
            # Check if we have data
            if self.module.get_data() is None:
                QMessageBox.warning(self, "Warning", "No data to export")
                return

            # Get file path
            file_path, _ = QFileDialog.getSaveFileName(
                self, "Export Plot", "", "PNG Files (*.png);;PDF Files (*.pdf);;SVG Files (*.svg)"
            )

            if not file_path:
                return

            # Save the figure
            self.canvas.fig.savefig(file_path, dpi=300, bbox_inches='tight')

            QMessageBox.information(self, "Export Successful", f"Plot exported to {file_path}")
            self.statusBar().showMessage(f"Plot exported to {file_path}", 2000)

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error exporting plot: {str(e)}")

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
        return self.module 
