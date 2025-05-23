#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=no-name-in-module, import-error
import sys
from pathlib import Path

# Add the root folder to the Python path
root_folder = Path(__file__).parent.parent.parent
sys.path.append(str(root_folder))

from PyQt5.QtWidgets import (
    QWidget,
    QGridLayout,
    QFormLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QDoubleSpinBox,
    QGroupBox,
    QMessageBox,
    QVBoxLayout,
    QFileDialog,
)
from PyQt5.QtCore import Qt, pyqtSlot, QMimeData
from PyQt5.QtGui import QDragEnterEvent, QDropEvent
import numpy as np
from packages.visualizer.coordinate_visualizer import CoordinateVisualizer
from packages.helpers import UnitConverter

class DragDropLineEdit(QLineEdit):
    """Custom QLineEdit that accepts drag and drop events."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setAcceptDrops(True)
        self.setPlaceholderText("Drag and drop CIF file here or click Browse...")
        self.setReadOnly(True)

    def dragEnterEvent(self, event: QDragEnterEvent):
        if event.mimeData().hasUrls():
            urls = event.mimeData().urls()
            if urls and urls[0].toLocalFile().endswith(".cif"):
                event.acceptProposedAction()

    def dropEvent(self, event: QDropEvent):
        urls = event.mimeData().urls()
        if urls:
            file_path = urls[0].toLocalFile()
            if file_path.endswith(".cif"):
                self.setText(file_path)
                self.textChanged.emit(file_path)


class InitWindow(QWidget):
    """Initialization window for setting up lattice parameters."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.parent = parent
        self.unit_converter = UnitConverter()
        self.init_ui()

    def init_ui(self):
        """Initialize UI components."""
        layout = QGridLayout(self)
        layout.setContentsMargins(20, 20, 20, 20)
        layout.setSpacing(40)

        # Group box for lattice parameters
        lattice_group = QGroupBox("Lattice Parameters")
        lattice_layout = QGridLayout(lattice_group)

        # Lattice constants (left column)
        self.a_input = QDoubleSpinBox()
        self.a_input.setRange(0.1, 100.0)
        self.a_input.setValue(5.0)
        self.a_input.setSuffix(" Å")
        lattice_layout.addWidget(QLabel("a:"), 0, 0)
        lattice_layout.addWidget(self.a_input, 0, 1)

        self.b_input = QDoubleSpinBox()
        self.b_input.setRange(0.1, 100.0)
        self.b_input.setValue(5.0)
        self.b_input.setSuffix(" Å")
        lattice_layout.addWidget(QLabel("b:"), 1, 0)
        lattice_layout.addWidget(self.b_input, 1, 1)

        self.c_input = QDoubleSpinBox()
        self.c_input.setRange(0.1, 100.0)
        self.c_input.setValue(5.0)
        self.c_input.setSuffix(" Å")
        lattice_layout.addWidget(QLabel("c:"), 2, 0)
        lattice_layout.addWidget(self.c_input, 2, 1)

        # Lattice angles (right column)
        self.alpha_input = QDoubleSpinBox()
        self.alpha_input.setRange(1.0, 179.0)
        self.alpha_input.setValue(90.0)
        self.alpha_input.setSuffix(" °")
        self.alpha_input.valueChanged.connect(self.update_visualization)
        lattice_layout.addWidget(QLabel("α:"), 0, 2)
        lattice_layout.addWidget(self.alpha_input, 0, 3)

        self.beta_input = QDoubleSpinBox()
        self.beta_input.setRange(1.0, 179.0)
        self.beta_input.setValue(90.0)
        self.beta_input.setSuffix(" °")
        self.beta_input.valueChanged.connect(self.update_visualization)
        lattice_layout.addWidget(QLabel("β:"), 1, 2)
        lattice_layout.addWidget(self.beta_input, 1, 3)

        self.gamma_input = QDoubleSpinBox()
        self.gamma_input.setRange(1.0, 179.0)
        self.gamma_input.setValue(90.0)
        self.gamma_input.setSuffix(" °")
        self.gamma_input.valueChanged.connect(self.update_visualization)
        lattice_layout.addWidget(QLabel("γ:"), 2, 2)
        lattice_layout.addWidget(self.gamma_input, 2, 3)

        # Add spacing between columns and margins
        lattice_layout.setColumnStretch(1, 1)
        lattice_layout.setColumnStretch(3, 1)
        lattice_layout.setHorizontalSpacing(40)
        lattice_layout.setVerticalSpacing(10)
        lattice_layout.setContentsMargins(20, 20, 20, 20)

        # Add lattice group to main layout at (0,0)
        layout.addWidget(lattice_group, 0, 0)

        # Group box for X-ray energy
        energy_group = QGroupBox("X-ray Energy")
        energy_layout = QFormLayout(energy_group)

        self.energy_input = QDoubleSpinBox()
        self.energy_input.setRange(1.0, 20000.0)
        self.energy_input.setValue(950.0)
        self.energy_input.setSuffix(" eV")
        self.energy_input.valueChanged.connect(self.on_energy_changed)
        energy_layout.addRow("Energy:", self.energy_input)

        self.wavelength_input = QDoubleSpinBox()
        self.wavelength_input.setRange(0.1, 100.0)
        self.wavelength_input.setValue(self.unit_converter.ev_to_angstrom(950.0))
        self.wavelength_input.setSuffix(" Å")
        self.wavelength_input.valueChanged.connect(self.on_wavelength_changed)
        energy_layout.addRow("Wavelength:", self.wavelength_input)

        # Add energy group to main layout at (0,1)
        layout.addWidget(energy_group, 0, 1)

        # Group box for Euler angles
        euler_group = QGroupBox("Euler Angles")
        euler_layout = QFormLayout(euler_group)

        self.roll_input = QDoubleSpinBox()
        self.roll_input.setRange(-180.0, 180.0)
        self.roll_input.setValue(0.0)
        self.roll_input.setSuffix(" °")
        self.roll_input.setToolTip("Rotation about the new X axis")
        self.roll_input.valueChanged.connect(self.update_visualization)
        euler_layout.addRow("Roll:", self.roll_input)

        self.pitch_input = QDoubleSpinBox()
        self.pitch_input.setRange(-180.0, 180.0)
        self.pitch_input.setValue(0.0)
        self.pitch_input.setSuffix(" °")
        self.pitch_input.setToolTip("Rotation about the new Y axis")
        self.pitch_input.valueChanged.connect(self.update_visualization)
        euler_layout.addRow("Pitch:", self.pitch_input)

        self.yaw_input = QDoubleSpinBox()
        self.yaw_input.setRange(-180.0, 180.0)
        self.yaw_input.setValue(0.0)
        self.yaw_input.setSuffix(" °")
        self.yaw_input.setToolTip("Rotation about the original Z axis")
        self.yaw_input.valueChanged.connect(self.update_visualization)
        euler_layout.addRow("Yaw:", self.yaw_input)

        # Add euler group to main layout at (0,2)
        layout.addWidget(euler_group, 0, 2)

        # Create and add the coordinate visualizer
        self.visualizer = CoordinateVisualizer()
        # initialize the visualizer
        self.visualizer.initialize(
            {
                "a": self.a_input.value(),
                "b": self.b_input.value(),
                "c": self.c_input.value(),
                "alpha": self.alpha_input.value(),
                "beta": self.beta_input.value(),
                "gamma": self.gamma_input.value(),
                "roll": self.roll_input.value(),
                "pitch": self.pitch_input.value(),
                "yaw": self.yaw_input.value(),
            }
        )
        self.visualizer.visualize_lab_system()
        layout.addWidget(self.visualizer, 1, 2)

        # File input area
        file_group = QGroupBox("Crystal Structure File")
        file_layout = QGridLayout(file_group)

        self.file_path_input = DragDropLineEdit()
        file_layout.addWidget(self.file_path_input, 0, 0)

        browse_button = QPushButton("Browse...")
        browse_button.clicked.connect(self.browse_cif_file)
        file_layout.addWidget(browse_button, 0, 1)

        # Add file group to main layout at (1,0) spanning 2 columns
        layout.addWidget(file_group, 1, 0, 1, 2)

        # Initialize button
        initialize_button = QPushButton("Initialize")
        initialize_button.clicked.connect(self.initialize)
        # Add initialize button at (2,0) spanning 2 columns
        layout.addWidget(initialize_button, 2, 0, 1, 2)

        # Add spacer
        layout.setRowStretch(3, 1)

    @pyqtSlot()
    def browse_cif_file(self):
        """Browse for CIF file."""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Open CIF File", "", "CIF Files (*.cif);;All Files (*)"
        )

        if file_path:
            self.file_path_input.setText(file_path)

    @pyqtSlot()
    def update_visualization(self):
        """Update the coordinate visualization when vectors change."""
        try:
            # Get current values
            roll = self.roll_input.value()
            pitch = self.pitch_input.value()
            yaw = self.yaw_input.value()

            # Validate values are within range
            if not (
                -180 <= roll <= 180 and -180 <= pitch <= 180 and -180 <= yaw <= 180
            ):
                # Reset to default values if invalid
                self.roll_input.setValue(0.0)
                self.pitch_input.setValue(0.0)
                self.yaw_input.setValue(0.0)
                roll, pitch, yaw = 0.0, 0.0, 0.0

            # Update the visualizer
            self.visualizer.initialize(
                {
                    "a": self.a_input.value(),
                    "b": self.b_input.value(),
                    "c": self.c_input.value(),
                    "alpha": self.alpha_input.value(),
                    "beta": self.beta_input.value(),
                    "gamma": self.gamma_input.value(),
                    "roll": roll,
                    "pitch": pitch,
                    "yaw": yaw,
                }
            )
            self.visualizer.visualize_lab_system()
        except Exception as e:
            raise e

    @pyqtSlot()
    def initialize(self):
        """Initialize the application with the provided parameters."""
        try:
            # Get parameters
            params = {
                "a": self.a_input.value(),
                "b": self.b_input.value(),
                "c": self.c_input.value(),
                "alpha": self.alpha_input.value(),
                "beta": self.beta_input.value(),
                "gamma": self.gamma_input.value(),
                "energy": self.energy_input.value(),
                "cif_file": (
                    self.file_path_input.text() if self.file_path_input.text() else None
                ),
                "roll": self.roll_input.value(),
                "pitch": self.pitch_input.value(),
                "yaw": self.yaw_input.value(),
            }

            # Pass parameters to parent (MainWindow)
            if self.parent:
                self.parent.set_parameters(params)
                self.parent.show_main_tabs()
                self.hide()

        except Exception as e:
            QMessageBox.critical(
                self, "Error", f"Error initializing parameters: {str(e)}"
            )

    @pyqtSlot()
    def on_energy_changed(self):
        """Update wavelength when energy changes."""
        try:
            # Block signals to prevent infinite loop
            self.wavelength_input.blockSignals(True)
            # Convert energy to wavelength
            wavelength = self.unit_converter.ev_to_angstrom(self.energy_input.value())
            self.wavelength_input.setValue(wavelength)
            # Unblock signals
            self.wavelength_input.blockSignals(False)
        except Exception as e:
            QMessageBox.warning(self, "Warning", f"Error converting energy: {str(e)}")

    @pyqtSlot()
    def on_wavelength_changed(self):
        """Update energy when wavelength changes."""
        try:
            # Block signals to prevent infinite loop
            self.energy_input.blockSignals(True)
            # Convert wavelength to energy
            energy = self.unit_converter.angstrom_to_ev(self.wavelength_input.value())
            self.energy_input.setValue(energy)
            # Unblock signals
            self.energy_input.blockSignals(False)
        except Exception as e:
            QMessageBox.warning(
                self, "Warning", f"Error converting wavelength: {str(e)}"
            )
