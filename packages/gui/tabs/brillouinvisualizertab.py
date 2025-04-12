#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
    QListWidgetItem,
    QCheckBox,
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
from packages.brillouin_visualizer.interface import BrillouinVisualizer


class MatplotlibCanvas(FigureCanvas):
    """Matplotlib canvas for embedding in Qt applications."""
    
    def __init__(self, width=6, height=5, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111, projection='3d')
        super(MatplotlibCanvas, self).__init__(self.fig)
        self.fig.tight_layout()


class BrillouinVisualizerTab(TabInterface):
    """Tab for visualizing Brillouin zones."""

    def __init__(self):
        # Create backend instance
        self.visualizer = BrillouinVisualizer()

        # Initialize UI
        super().__init__()

        # Set window title
        self.setWindowTitle("Brillouin Zone Visualizer")

    def init_ui(self):
        """Initialize UI components."""
        # Create main layout with two columns
        main_layout = QGridLayout()
        self.layout.addLayout(main_layout)

        # Left panel - controls
        left_panel = QGridLayout()
        main_layout.addLayout(left_panel, 0, 0)

        # Right panel - visualization
        right_panel = QGridLayout()
        main_layout.addLayout(right_panel, 0, 1)

        # Create left panel components
        self.create_structure_input(left_panel)
        self.create_view_controls(left_panel)
        self.create_path_controls(left_panel)

        # Create visualization area
        self.create_visualization_area(right_panel)

    def create_structure_input(self, parent_layout):
        """Create input area for structure parameters."""
        structure_group = QGroupBox("Crystal Structure")
        structure_layout = QGridLayout(structure_group)

        # Structure input tabs
        input_tabs = QTabWidget()
        structure_layout.addWidget(input_tabs, 0, 0)

        # File input tab
        file_tab = QWidget()
        file_layout = QGridLayout(file_tab)

        file_form = QFormLayout()
        self.file_path_input = QLineEdit()
        self.file_path_input.setPlaceholderText("No file selected")
        self.file_path_input.setReadOnly(True)
        browse_button = QPushButton("Browse...")
        browse_button.clicked.connect(self.browse_cif_file)

        file_form.addRow("CIF File:", self.file_path_input)
        file_layout.addLayout(file_form, 0, 0)
        file_layout.addWidget(browse_button, 0, 1)

        load_file_button = QPushButton("Load Structure from File")
        load_file_button.clicked.connect(self.load_from_file)
        file_layout.addWidget(load_file_button, 1, 0, 1, 2)

        input_tabs.addTab(file_tab, "From File")

        # Manual input tab
        manual_tab = QWidget()
        manual_layout = QGridLayout(manual_tab)

        # Lattice type selection
        lattice_form = QFormLayout()
        self.lattice_type_combo = QComboBox()
        self.lattice_type_combo.addItems(["Simple Cubic (SC)", "Body-Centered Cubic (BCC)", 
                                        "Face-Centered Cubic (FCC)", "Hexagonal", "Custom"])
        lattice_form.addRow("Lattice Type:", self.lattice_type_combo)
        manual_layout.addLayout(lattice_form, 0, 0)

        # Lattice constants group
        constants_group = QGroupBox("Lattice Constants")
        constants_layout = QFormLayout(constants_group)

        self.a_input = QDoubleSpinBox()
        self.a_input.setRange(0.1, 100.0)
        self.a_input.setValue(5.0)
        self.a_input.setSuffix(" Å")
        constants_layout.addRow("a:", self.a_input)

        self.b_input = QDoubleSpinBox()
        self.b_input.setRange(0.1, 100.0)
        self.b_input.setValue(5.0)
        self.b_input.setSuffix(" Å")
        constants_layout.addRow("b:", self.b_input)

        self.c_input = QDoubleSpinBox()
        self.c_input.setRange(0.1, 100.0)
        self.c_input.setValue(5.0)
        self.c_input.setSuffix(" Å")
        constants_layout.addRow("c:", self.c_input)

        self.alpha_input = QDoubleSpinBox()
        self.alpha_input.setRange(1.0, 179.0)
        self.alpha_input.setValue(90.0)
        self.alpha_input.setSuffix(" °")
        constants_layout.addRow("α:", self.alpha_input)

        self.beta_input = QDoubleSpinBox()
        self.beta_input.setRange(1.0, 179.0)
        self.beta_input.setValue(90.0)
        self.beta_input.setSuffix(" °")
        constants_layout.addRow("β:", self.beta_input)

        self.gamma_input = QDoubleSpinBox()
        self.gamma_input.setRange(1.0, 179.0)
        self.gamma_input.setValue(90.0)
        self.gamma_input.setSuffix(" °")
        constants_layout.addRow("γ:", self.gamma_input)

        manual_layout.addWidget(constants_group, 1, 0)

        # Connect lattice type changes to enable/disable fields
        self.lattice_type_combo.currentIndexChanged.connect(self.update_lattice_inputs)

        # Load button
        load_manual_button = QPushButton("Generate Brillouin Zone")
        load_manual_button.clicked.connect(self.load_from_parameters)
        manual_layout.addWidget(load_manual_button, 2, 0)

        input_tabs.addTab(manual_tab, "Manual Input")

        # Add to parent layout
        parent_layout.addWidget(structure_group, 0, 0)

        # Initial update of input fields
        self.update_lattice_inputs(0)

    def create_view_controls(self, parent_layout):
        """Create controls for view options."""
        view_group = QGroupBox("View Options")
        view_layout = QGridLayout(view_group)

        # View direction selection
        direction_form = QFormLayout()
        self.view_direction_combo = QComboBox()
        self.view_direction_combo.addItems(["3D View", "View along x", "View along y", "View along z"])
        self.view_direction_combo.currentIndexChanged.connect(self.update_visualization)
        direction_form.addRow("View Direction:", self.view_direction_combo)
        view_layout.addLayout(direction_form, 0, 0)

        # Special points checkbox
        self.show_points_check = QCheckBox("Show Special Points")
        self.show_points_check.setChecked(True)
        self.show_points_check.stateChanged.connect(self.update_visualization)
        view_layout.addWidget(self.show_points_check, 1, 0)

        parent_layout.addWidget(view_group, 1, 0)

    def create_path_controls(self, parent_layout):
        """Create controls for path visualization."""
        path_group = QGroupBox("Path in Brillouin Zone")
        path_layout = QGridLayout(path_group)

        # Special points list
        points_label = QLabel("Special Points:")
        path_layout.addWidget(points_label, 0, 0)

        self.points_list = QListWidget()
        self.points_list.setSelectionMode(QListWidget.ExtendedSelection)
        path_layout.addWidget(self.points_list)

        # Path controls
        path_buttons_layout = QGridLayout()

        add_point_button = QPushButton("Add to Path")
        add_point_button.clicked.connect(self.add_to_path)
        path_buttons_layout.addWidget(add_point_button, 0, 0)

        clear_path_button = QPushButton("Clear Path")
        clear_path_button.clicked.connect(self.clear_path)
        path_buttons_layout.addWidget(clear_path_button, 0, 1)

        path_layout.addLayout(path_buttons_layout, 1, 0)

        # Current path
        path_layout.addWidget(QLabel("Current Path:"), 2, 0)
        self.path_list = QListWidget()
        path_layout.addWidget(self.path_list)

        # Visualize path button
        visualize_path_button = QPushButton("Visualize Path")
        visualize_path_button.clicked.connect(self.visualize_path)
        path_layout.addWidget(visualize_path_button)

        parent_layout.addWidget(path_group)

    def create_visualization_area(self, parent_layout):
        """Create area for visualization."""
        vis_group = QGroupBox("Visualization")
        vis_layout = QGridLayout(vis_group)

        # Add matplotlib canvas
        self.canvas = MatplotlibCanvas(width=8, height=6)
        vis_layout.addWidget(self.canvas, 0, 0)

        parent_layout.addWidget(vis_group)

    @pyqtSlot()
    def browse_cif_file(self):
        """Browse for CIF file."""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Open CIF File", "", "CIF Files (*.cif);;All Files (*)"
        )

        if file_path:
            self.file_path_input.setText(file_path)

    @pyqtSlot()
    def load_from_file(self):
        """Load structure from CIF file."""
        file_path = self.file_path_input.text()
        if not file_path:
            QMessageBox.warning(self, "Warning", "Please select a CIF file first.")
            return

        try:
            success = self.visualizer.initialize_from_cif(file_path)

            if success:
                QMessageBox.information(self, "Success", "Crystal structure loaded successfully!")
                self.update_structure_info()
                self.update_special_points_list()
                self.update_visualization()
            else:
                QMessageBox.warning(self, "Error", "Failed to load crystal structure.")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error loading crystal structure: {str(e)}")

    @pyqtSlot()
    def load_from_parameters(self):
        """Load structure from manual parameters."""
        try:
            # Get lattice type
            lattice_type_text = self.lattice_type_combo.currentText()
            if "SC" in lattice_type_text:
                lattice_type = "SC"
            elif "BCC" in lattice_type_text:
                lattice_type = "BCC"
            elif "FCC" in lattice_type_text:
                lattice_type = "FCC"
            elif "Hexagonal" in lattice_type_text:
                lattice_type = "HEX"
            else:
                lattice_type = "CUSTOM"

            # Get lattice constants
            lattice_constants = {
                'a': self.a_input.value(),
                'b': self.b_input.value(),
                'c': self.c_input.value(),
                'alpha': self.alpha_input.value(),
                'beta': self.beta_input.value(),
                'gamma': self.gamma_input.value()
            }

            success = self.visualizer.initialize_from_parameters(lattice_type, lattice_constants)

            if success:
                QMessageBox.information(self, "Success", "Brillouin zone generated successfully!")
                self.update_structure_info()
                self.update_special_points_list()
                self.update_visualization()
            else:
                QMessageBox.warning(self, "Error", "Failed to generate Brillouin zone.")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error generating Brillouin zone: {str(e)}")

    @pyqtSlot(int)
    def update_lattice_inputs(self, index):
        """Update lattice constant inputs based on selected lattice type."""
        is_custom = index == 4  # "Custom" is the fifth option (index 4)
        is_hexagonal = index == 3  # "Hexagonal" is the fourth option (index 3)

        # For non-custom lattices, some parameters are fixed
        if not is_custom:
            if is_hexagonal:
                # Hexagonal: a=b≠c, alpha=beta=90°, gamma=120°
                self.b_input.setValue(self.a_input.value())
                self.b_input.setEnabled(False)
                self.alpha_input.setValue(90.0)
                self.alpha_input.setEnabled(False)
                self.beta_input.setValue(90.0)
                self.beta_input.setEnabled(False)
                self.gamma_input.setValue(120.0)
                self.gamma_input.setEnabled(False)
            else:
                # Cubic lattices: a=b=c, all angles 90°
                self.b_input.setValue(self.a_input.value())
                self.b_input.setEnabled(False)
                self.c_input.setValue(self.a_input.value())
                self.c_input.setEnabled(False)
                self.alpha_input.setValue(90.0)
                self.alpha_input.setEnabled(False)
                self.beta_input.setValue(90.0)
                self.beta_input.setEnabled(False)
                self.gamma_input.setValue(90.0)
                self.gamma_input.setEnabled(False)
        else:
            # Custom: all parameters editable
            self.b_input.setEnabled(True)
            self.c_input.setEnabled(True)
            self.alpha_input.setEnabled(True)
            self.beta_input.setEnabled(True)
            self.gamma_input.setEnabled(True)

    def update_structure_info(self):
        """Update UI with information about the loaded structure."""
        structure = self.visualizer.get_crystal_structure()
        if structure:
            # In a real implementation, this might update labels or other UI elements
            # with information about the crystal structure
            pass

    def update_special_points_list(self):
        """Update the list of special points."""
        self.points_list.clear()
        special_points = self.visualizer.get_special_points()

        for name, coords in special_points.items():
            item_text = f"{name}: ({coords[0]:.3f}, {coords[1]:.3f}, {coords[2]:.3f})"
            item = QListWidgetItem(item_text)
            item.setData(Qt.UserRole, name)  # Store point name for reference
            self.points_list.addItem(item)

    @pyqtSlot()
    def update_visualization(self):
        """Update the visualization based on current settings."""
        if not self.visualizer.is_initialized():
            return

        try:
            # Get view direction
            view_index = self.view_direction_combo.currentIndex()
            if view_index == 0:
                view_direction = None  # 3D view
            elif view_index == 1:
                view_direction = 'x'
            elif view_index == 2:
                view_direction = 'y'
            else:
                view_direction = 'z'

            # Get whether to show special points
            show_points = self.show_points_check.isChecked()

            # Generate figure
            fig = self.visualizer.visualize_brillouin_zone(show_points, view_direction)

            if fig:
                # Replace the current figure in the canvas
                self.canvas.fig.clear()

                # If current view is 3D but new view is 2D or vice versa, need to recreate axes
                old_projection = getattr(self.canvas.axes, 'projection', None)
                new_projection = '3d' if view_direction is None else None

                if old_projection != new_projection:
                    # Remove old axes and create new ones with correct projection
                    self.canvas.fig.delaxes(self.canvas.axes)
                    if new_projection:
                        self.canvas.axes = self.canvas.fig.add_subplot(111, projection=new_projection)
                    else:
                        self.canvas.axes = self.canvas.fig.add_subplot(111)

                # Copy content from generated figure to canvas figure
                for i, ax in enumerate(fig.get_axes()):
                    # Get content from the generated figure axes
                    if i == 0:  # Only handling the first subplot for simplicity
                        # Copy properties to existing axes
                        if new_projection == '3d':
                            # 3D axes
                            for collection in ax.collections:
                                self.canvas.axes.add_collection3d(collection.copy())
                            for line in ax.lines:
                                self.canvas.axes.plot3D(line.get_xdata(), line.get_ydata(), line.get_zdata(),
                                                     color=line.get_color(), linestyle=line.get_linestyle(),
                                                     linewidth=line.get_linewidth(), alpha=line.get_alpha())
                            for text in ax.texts:
                                self.canvas.axes.text(text.get_position()[0], text.get_position()[1], 
                                                   text.get_position()[2], text.get_text(),
                                                   fontsize=text.get_fontsize(), color=text.get_color())

                            # Copy axis settings
                            self.canvas.axes.set_xlim(ax.get_xlim())
                            self.canvas.axes.set_ylim(ax.get_ylim())
                            self.canvas.axes.set_zlim(ax.get_zlim())
                            self.canvas.axes.set_xlabel(ax.get_xlabel())
                            self.canvas.axes.set_ylabel(ax.get_ylabel())
                            self.canvas.axes.set_zlabel(ax.get_zlabel())
                            self.canvas.axes.set_title(ax.get_title())
                        else:
                            # 2D axes
                            for collection in ax.collections:
                                self.canvas.axes.add_collection(collection.copy())
                            for line in ax.lines:
                                self.canvas.axes.plot(line.get_xdata(), line.get_ydata(),
                                                   color=line.get_color(), linestyle=line.get_linestyle(),
                                                   linewidth=line.get_linewidth(), alpha=line.get_alpha())
                            for text in ax.texts:
                                self.canvas.axes.text(text.get_position()[0], text.get_position()[1], 
                                                   text.get_text(), fontsize=text.get_fontsize(), 
                                                   color=text.get_color())

                            # Copy axis settings
                            self.canvas.axes.set_xlim(ax.get_xlim())
                            self.canvas.axes.set_ylim(ax.get_ylim())
                            self.canvas.axes.set_xlabel(ax.get_xlabel())
                            self.canvas.axes.set_ylabel(ax.get_ylabel())
                            self.canvas.axes.set_title(ax.get_title())
                            if hasattr(ax, 'get_aspect') and ax.get_aspect() == 'equal':
                                self.canvas.axes.set_aspect('equal')

                self.canvas.draw()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error updating visualization: {str(e)}")

    @pyqtSlot()
    def add_to_path(self):
        """Add selected special points to the path."""
        selected_items = self.points_list.selectedItems()
        if not selected_items:
            return

        for item in selected_items:
            point_name = item.data(Qt.UserRole)
            path_item = QListWidgetItem(point_name)
            path_item.setData(Qt.UserRole, point_name)
            self.path_list.addItem(path_item)

    @pyqtSlot()
    def clear_path(self):
        """Clear the current path."""
        self.path_list.clear()

    @pyqtSlot()
    def visualize_path(self):
        """Visualize the selected path through the Brillouin zone."""
        if not self.visualizer.is_initialized():
            QMessageBox.warning(self, "Warning", "Please load a crystal structure first.")
            return

        if self.path_list.count() < 2:
            QMessageBox.warning(self, "Warning", "Please select at least two points for the path.")
            return

        try:
            # Get path points
            path = []
            for i in range(self.path_list.count()):
                item = self.path_list.item(i)
                point_name = item.data(Qt.UserRole)
                path.append(point_name)

            # Generate path visualization
            path_fig, _ = self.visualizer.visualize_path(path)

            if path_fig:
                # Replace the current figure in the canvas
                self.canvas.fig.clear()

                # Ensure we have 3D axes
                if not hasattr(self.canvas.axes, 'projection') or self.canvas.axes.projection != '3d':
                    self.canvas.fig.delaxes(self.canvas.axes)
                    self.canvas.axes = self.canvas.fig.add_subplot(111, projection='3d')

                # Copy content from generated figure to canvas figure
                for i, ax in enumerate(path_fig.get_axes()):
                    if i == 0:  # Only handling the first subplot for simplicity
                        # Copy 3D content
                        for collection in ax.collections:
                            self.canvas.axes.add_collection3d(collection.copy())
                        for line in ax.lines:
                            self.canvas.axes.plot3D(line.get_xdata(), line.get_ydata(), line.get_zdata(),
                                                 color=line.get_color(), linestyle=line.get_linestyle(),
                                                 linewidth=line.get_linewidth(), alpha=line.get_alpha())
                        for text in ax.texts:
                            self.canvas.axes.text(text.get_position()[0], text.get_position()[1], 
                                               text.get_position()[2], text.get_text(),
                                               fontsize=text.get_fontsize(), color=text.get_color())

                        # Copy axis settings
                        self.canvas.axes.set_xlim(ax.get_xlim())
                        self.canvas.axes.set_ylim(ax.get_ylim())
                        self.canvas.axes.set_zlim(ax.get_zlim())
                        self.canvas.axes.set_xlabel(ax.get_xlabel())
                        self.canvas.axes.set_ylabel(ax.get_ylabel())
                        self.canvas.axes.set_zlabel(ax.get_zlabel())
                        self.canvas.axes.set_title(ax.get_title())

                self.canvas.draw()

                # Update view direction combo to 3D
                self.view_direction_combo.setCurrentIndex(0)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error visualizing path: {str(e)}")

    def open_file(self, file_path):
        """Handle opening a file from the main window."""
        if file_path.lower().endswith('.cif'):
            self.file_path_input.setText(file_path)
            self.load_from_file()
            return True
        return False

    def get_module_instance(self):
        """Get the backend module instance."""
        return self.visualizer 
