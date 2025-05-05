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
        lattice_layout.addWidget(QLabel("α:"), 0, 2)
        lattice_layout.addWidget(self.alpha_input, 0, 3)

        self.beta_input = QDoubleSpinBox()
        self.beta_input.setRange(1.0, 179.0)
        self.beta_input.setValue(90.0)
        self.beta_input.setSuffix(" °")
        lattice_layout.addWidget(QLabel("β:"), 1, 2)
        lattice_layout.addWidget(self.beta_input, 1, 3)

        self.gamma_input = QDoubleSpinBox()
        self.gamma_input.setRange(1.0, 179.0)
        self.gamma_input.setValue(90.0)
        self.gamma_input.setSuffix(" °")
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
        energy_layout.addRow("Energy:", self.energy_input)

        # Add energy group to main layout at (0,1)
        layout.addWidget(energy_group, 0, 1)

        # Group box for coordinate system
        coord_group = QGroupBox("Coordinate System")
        coord_layout = QFormLayout(coord_group)

        self.e_H_input = QLineEdit()
        self.e_H_input.setPlaceholderText("1 0 0")
        self.e_H_input.setToolTip(
            "Enter three numbers separated by spaces (e.g. 0 -1 1)"
        )
        self.e_H_input.textChanged.connect(self.update_visualization)
        coord_layout.addRow("e_H:", self.e_H_input)

        self.e_K_input = QLineEdit()
        self.e_K_input.setPlaceholderText("0 1 0")
        self.e_K_input.setToolTip(
            "Enter three numbers separated by spaces (e.g. 0 -1 1)"
        )
        self.e_K_input.textChanged.connect(self.update_visualization)
        coord_layout.addRow("e_K:", self.e_K_input)

        self.e_L_input = QLineEdit()
        self.e_L_input.setPlaceholderText("0 0 1")
        self.e_L_input.setToolTip(
            "Enter three numbers separated by spaces (e.g. 0 -1 1)"
        )
        self.e_L_input.textChanged.connect(self.update_visualization)
        coord_layout.addRow("e_L:", self.e_L_input)

        # Add coord group to main layout at (0,2)
        layout.addWidget(coord_group, 0, 2)

        # Create and add the coordinate visualizer
        self.visualizer = CoordinateVisualizer()
        # initialize the visualizer
        self.visualizer.visualize_lab_system(
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
        )
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
            e_H = parse_vector(self.e_H_input.text(), default=[1, 0, 0])
            e_K = parse_vector(self.e_K_input.text(), default=[0, 1, 0])
            e_L = parse_vector(self.e_L_input.text(), default=[0, 0, 1])

            # Update the visualization
            self.visualizer.visualize_lab_system(e_H, e_K, e_L)
        except Exception as e:
            # If there's an error in parsing, just keep the current visualization
            pass

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
                "e_H": parse_vector(self.e_H_input.text(), default=[1, 0, 0]),
                "e_K": parse_vector(self.e_K_input.text(), default=[0, 1, 0]),
                "e_L": parse_vector(self.e_L_input.text(), default=[0, 0, 1]),
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


def parse_vector(input_str, default=[0, 0, 1]):
    try:
        # Split on whitespace and convert to float and normalize it
        components = [float(x) for x in input_str.strip().split()]
        if len(components) != 3:
            return np.array(default)
        return np.array(components) / np.linalg.norm(components)
    except:
        return np.array(default)
