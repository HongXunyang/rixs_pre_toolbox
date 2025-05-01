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
    QVBoxLayout,
    QHBoxLayout,
)
from PyQt5.QtCore import Qt, pyqtSlot, QMimeData
from PyQt5.QtGui import QDragEnterEvent, QDropEvent
import sys
import os
import matplotlib

# Tell matplotlib to render plots using the Qt5 framework with the Agg (Anti-Grain Geometry) backend for drawing
matplotlib.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

# Add parent directory to path to allow imports from sibling packages
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from packages.gui.tabs.tab_interface import TabInterface
from packages.brillouin_calculator.interface import BrillouinCalculator
from packages.visualizer.scattering_visualizer import ScatteringVisualizer
from packages.helpers.tips import Tips, set_tip

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
                # Emit the textChanged signal to notify parent
                self.textChanged.emit(file_path)


class BrillouinCalculatorTab(TabInterface):
    """Tab for calculating Brillouin zone parameters."""

    def __init__(self, main_window=None):
        # Create backend instance
        self.calculator = BrillouinCalculator()
        self.visualizer = ScatteringVisualizer()
        self.hkl_visualizer = ScatteringVisualizer()
        self.tips = Tips()

        # Initialize UI
        super().__init__(main_window)
        self.init_ui()
        params = self.main_window.get_lattice_parameters()
        self.set_lattice_parameters(params)
        # Set window title
        self.setWindowTitle("Brillouin Zone Calculator")

    def init_ui(self):
        """Initialize UI components."""
        # Create tab widget for input methods
        self.tab_widget = QTabWidget()
        self.layout.addWidget(self.tab_widget, 0, 0)  # Add to grid at (0,0)

        # Create tabs for different functionalities
        self.create_angles_to_hkl_tab()
        self.create_hkl_to_angles_tab()

    def set_tip(self, widget, name):
        """Set the tooltip and status tip for a widget by the name"""
        set_tip(widget, self.tips.tip(name))

    def set_lattice_parameters(self, params):
        """Set lattice parameters from global settings."""
        try:

            success_1 = self.calculator.initialize(
                a=params["a"],
                b=params["b"],
                c=params["c"],
                alpha=params["alpha"],
                beta=params["beta"],
                gamma=params["gamma"],
                energy=params["energy"],
                e_H=params["e_H"],
                e_K=params["e_K"],
                e_L=params["e_L"],
            )

            self.visualizer.initialize(
                e_H=params["e_H"],
                e_K=params["e_K"],
                e_L=params["e_L"],
            )
            self.hkl_visualizer.initialize(
                e_H=params["e_H"],
                e_K=params["e_K"],
                e_L=params["e_L"],
            )
            if not success_1:
                QMessageBox.warning(self, "Error", "Failed to initialize calculator!")

        except Exception as e:
            QMessageBox.critical(
                self, "Error", f"Error initializing calculator: {str(e)}"
            )

    def create_angles_to_hkl_tab(self):
        """Create tab for angles to HKL calculation."""
        angles_tab = QWidget()
        angles_layout = QGridLayout(angles_tab)

        # Input form
        form_group = QGroupBox("Scattering Angles")
        form_layout = QFormLayout(form_group)

        self.tth_angle_input = QDoubleSpinBox()
        self.tth_angle_input.setRange(0.0, 180.0)
        self.tth_angle_input.setValue(150.0)
        self.tth_angle_input.setSuffix(" °")
        self.set_tip(self.tth_angle_input, "TTH")
        form_layout.addRow("tth:", self.tth_angle_input)

        self.theta_angle_input = QDoubleSpinBox()
        self.theta_angle_input.setRange(-180.0, 180.0)
        self.theta_angle_input.setValue(50.0)
        self.theta_angle_input.setSuffix(" °")
        self.set_tip(self.theta_angle_input, "THETA")
        form_layout.addRow("θ:", self.theta_angle_input)

        self.phi_angle_input = QDoubleSpinBox()
        self.phi_angle_input.setRange(-180.0, 180.0)
        self.phi_angle_input.setValue(0.0)
        self.phi_angle_input.setSuffix(" °")
        form_layout.addRow("φ:", self.phi_angle_input)

        self.chi_angle_input = QDoubleSpinBox()
        self.chi_angle_input.setRange(-180.0, 180.0)
        self.chi_angle_input.setValue(0.0)
        self.chi_angle_input.setSuffix(" °")
        form_layout.addRow("χ:", self.chi_angle_input)

        angles_layout.addWidget(form_group, 0, 0)

        # Calculate button
        calculate_button = QPushButton("Calculate HKL")
        calculate_button.clicked.connect(self.calculate_hkl)
        angles_layout.addWidget(calculate_button, 1, 0)

        # Results group
        results_group = QGroupBox("Results")
        results_layout = QFormLayout(results_group)

        self.H_result = QLineEdit()
        self.H_result.setReadOnly(True)
        results_layout.addRow("H:", self.H_result)

        self.K_result = QLineEdit()
        self.K_result.setReadOnly(True)
        results_layout.addRow("K:", self.K_result)

        self.L_result = QLineEdit()
        self.L_result.setReadOnly(True)
        results_layout.addRow("L:", self.L_result)

        angles_layout.addWidget(results_group, 2, 0)

        # Visualizer
        self.visualizer.visualize_lab_system(is_clear=True)
        self.visualizer.visualize_scattering_geometry(is_clear=False)
        angles_layout.addWidget(self.visualizer, 0, 1, 2, 1)

        # Add to tab widget
        self.tab_widget.addTab(angles_tab, "Angles → HKL")

    def create_hkl_to_angles_tab(self):
        """Create tab for HKL to angles calculation."""
        hkl_tab = QWidget()
        hkl_layout = QGridLayout(hkl_tab)

        # Input form
        form_group = QGroupBox("HKL Indices")
        form_layout = QFormLayout(form_group)

        self.H_input = QDoubleSpinBox()
        self.H_input.setRange(-10.0, 10.0)
        self.H_input.setDecimals(3)
        self.H_input.setValue(0.15)
        form_layout.addRow("H:", self.H_input)

        self.K_input = QDoubleSpinBox()
        self.K_input.setRange(-10.0, 10.0)
        self.K_input.setDecimals(3)
        self.K_input.setValue(0.1)
        form_layout.addRow("K:", self.K_input)

        self.L_input = QDoubleSpinBox()
        self.L_input.setRange(-10.0, 10.0)
        self.L_input.setDecimals(3)
        self.L_input.setValue(-0.5)
        form_layout.addRow("L:", self.L_input)

        hkl_layout.addWidget(form_group, 0, 0)

        # Constraints group
        constraints_group = QGroupBox("Constraints")
        constraints_layout = QFormLayout(constraints_group)

        self.tth_max_input = QDoubleSpinBox()
        self.tth_max_input.setRange(0.0, 180.0)
        self.tth_max_input.setValue(152.0)
        self.tth_max_input.setSuffix(" °")
        constraints_layout.addRow("tth max:", self.tth_max_input)

        # Fix angle selection
        fix_angle_group = QGroupBox("Fix Angle")
        fix_angle_layout = QHBoxLayout(fix_angle_group)

        self.fix_chi_radio = QRadioButton("Fix χ")
        self.fix_phi_radio = QRadioButton("Fix φ")
        self.fix_chi_radio.setChecked(True)  # Default to fixed chi

        fix_angle_layout.addWidget(self.fix_chi_radio)
        fix_angle_layout.addWidget(self.fix_phi_radio)

        constraints_layout.addRow(fix_angle_group)

        # constaint of chi
        self.chi_input = QDoubleSpinBox()
        self.chi_input.setRange(-180.0, 180.0)
        self.chi_input.setValue(0.0)
        self.chi_input.setSuffix(" °")
        constraints_layout.addRow("χ:", self.chi_input)

        # constraint of phi
        self.phi_input = QDoubleSpinBox()
        self.phi_input.setRange(-180.0, 180.0)
        self.phi_input.setValue(0.0)
        self.phi_input.setSuffix(" °")
        constraints_layout.addRow("φ:", self.phi_input)

        # Connect radio buttons to enable/disable corresponding inputs
        self.fix_chi_radio.toggled.connect(self.update_fixed_angle_ui)
        self.fix_phi_radio.toggled.connect(self.update_fixed_angle_ui)

        # Initialize UI state
        self.update_fixed_angle_ui()

        hkl_layout.addWidget(constraints_group, 1, 0)

        # Calculate button
        calculate_button = QPushButton("Calculate Angles")
        calculate_button.clicked.connect(self.calculate_angles)
        hkl_layout.addWidget(calculate_button, 2, 0)

        # Results group
        results_group = QGroupBox("Results")
        results_layout = QVBoxLayout(results_group)

        # First row: tth and phi
        first_row = QWidget()
        first_row_layout = QHBoxLayout(first_row)
        first_row_layout.setContentsMargins(0, 0, 0, 0)

        tth_widget = QWidget()
        tth_layout = QFormLayout(tth_widget)
        tth_layout.setContentsMargins(0, 0, 0, 0)
        self.tth_result = QLineEdit()
        self.set_tip(self.tth_result, "TTH")
        self.tth_result.setReadOnly(True)
        tth_layout.addRow("tth:", self.tth_result)
        first_row_layout.addWidget(tth_widget)

        phi_widget = QWidget()
        phi_layout = QFormLayout(phi_widget)
        phi_layout.setContentsMargins(0, 0, 0, 0)
        self.phi_result = QLineEdit()
        self.phi_result.setReadOnly(True)
        phi_layout.addRow("φ:", self.phi_result)
        first_row_layout.addWidget(phi_widget)

        results_layout.addWidget(first_row)

        # Second row: theta and chi
        second_row = QWidget()
        second_row_layout = QHBoxLayout(second_row)
        second_row_layout.setContentsMargins(0, 0, 0, 0)

        theta_widget = QWidget()
        theta_layout = QFormLayout(theta_widget)
        theta_layout.setContentsMargins(0, 0, 0, 0)
        self.theta_result = QLineEdit()
        self.set_tip(self.theta_result, "THETA")
        self.theta_result.setReadOnly(True)
        theta_layout.addRow("θ:", self.theta_result)
        second_row_layout.addWidget(theta_widget)

        chi_widget = QWidget()
        chi_layout = QFormLayout(chi_widget)
        chi_layout.setContentsMargins(0, 0, 0, 0)
        self.chi_result = QLineEdit()
        self.chi_result.setReadOnly(True)
        chi_layout.addRow("χ:", self.chi_result)
        second_row_layout.addWidget(chi_widget)

        results_layout.addWidget(second_row)

        hkl_layout.addWidget(results_group, 3, 0)

        # Visualizer
        self.hkl_visualizer.visualize_lab_system(is_clear=True)
        self.hkl_visualizer.visualize_scattering_geometry(is_clear=False)
        hkl_layout.addWidget(self.hkl_visualizer, 0, 1, 3, 1)

        # Add to tab widget
        self.tab_widget.addTab(hkl_tab, "HKL → Angles")

    @pyqtSlot()
    def browse_cif_file(self):
        """Browse for CIF file."""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Open CIF File", "", "CIF Files (*.cif);;All Files (*)"
        )

        if file_path:
            self.file_path_input.setText(file_path)

    @pyqtSlot()
    def calculate_hkl(self):
        """Calculate HKL from angles."""
        try:
            # Check if calculator is initialized
            if not self.calculator.is_initialized():
                QMessageBox.warning(
                    self, "Warning", "Please initialize the calculator first!"
                )
                self.tab_widget.setCurrentIndex(0)
                return

            # Calculate HKL
            tth = self.tth_angle_input.value()
            theta = self.theta_angle_input.value()
            phi = self.phi_angle_input.value()
            chi = self.chi_angle_input.value()
            result = self.calculator.calculate_hkl(
                tth=tth,
                theta=theta,
                phi=phi,
                chi=chi,
            )
            success = result.get("success", False)
            if not success:
                QMessageBox.warning(
                    self, "Warning", result.get("error", "Unknown error")
                )
                return
            # Update results
            self.H_result.setText(f"{result['H']:.4f}")
            self.K_result.setText(f"{result['K']:.4f}")
            self.L_result.setText(f"{result['L']:.4f}")

            # Update visualization
            self.visualizer.visualize_lab_system(is_clear=True, chi=chi, phi=phi)
            self.visualizer.visualize_scattering_geometry(
                scattering_angles=result, is_clear=False
            )

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error calculating HKL: {str(e)}")

    @pyqtSlot()
    def update_fixed_angle_ui(self):
        """Update UI based on which angle is fixed."""
        is_chi_fixed = self.fix_chi_radio.isChecked()
        self.chi_input.setEnabled(is_chi_fixed)
        self.phi_input.setEnabled(not is_chi_fixed)

    @pyqtSlot()
    def calculate_angles(self):
        """Calculate angles from HKL."""
        try:
            # Check if calculator is initialized
            if not self.calculator.is_initialized():
                QMessageBox.warning(
                    self, "Warning", "Please initialize the calculator first!"
                )
                self.tab_widget.setCurrentIndex(0)
                return

            # Determine which angle to fix based on radio button selection
            fixed_angle_name = "chi" if self.fix_chi_radio.isChecked() else "phi"
            fixed_angle_value = (
                self.chi_input.value()
                if self.fix_chi_radio.isChecked()
                else self.phi_input.value()
            )

            # Calculate angles
            result = self.calculator.calculate_angles(
                H_crystal=self.H_input.value(),
                K_crystal=self.K_input.value(),
                L_crystal=self.L_input.value(),
                fixed_angle=fixed_angle_value,
                fixed_angle_name=fixed_angle_name,
                tth_max=self.tth_max_input.value(),
            )

            success = result.get("success", False)
            if not success:
                QMessageBox.warning(
                    self, "Warning", result.get("error", "Unknown error")
                )
                return

            # Update results
            tth = result["tth"]
            theta = result["theta"]
            chi = result["chi"]
            phi = result["phi"]
            self.tth_result.setText(f"{tth:.4f}°")
            self.theta_result.setText(f"{theta:.4f}°")
            self.phi_result.setText(f"{phi:.4f}°")
            self.chi_result.setText(f"{chi:.4f}°")

            # Update visualization
            # self.update_visualization(scattering_angles=result)
            # Update hkl tab visualization
            self.hkl_visualizer.visualize_lab_system(is_clear=True, chi=chi, phi=phi)
            self.hkl_visualizer.visualize_scattering_geometry(
                scattering_angles=result, is_clear=False
            )

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error calculating angles: {str(e)}")

    def get_module_instance(self):
        """Get the backend module instance."""
        return self.calculator

    def get_state(self):
        """Get the current state for session saving."""
        return {
            "lattice": {
                "a": self.a_input.value(),
                "b": self.b_input.value(),
                "c": self.c_input.value(),
                "alpha": self.alpha_input.value(),
                "beta": self.beta_input.value(),
                "gamma": self.gamma_input.value(),
            },
            "energy": self.energy_input.value(),
            "file_path": self.file_path_input.text(),
            "current_tab": self.tab_widget.currentIndex(),
        }

    def set_state(self, state):
        """Restore tab state from saved session."""
        try:
            if "lattice" in state:
                lattice = state["lattice"]
                self.a_input.setValue(lattice.get("a", 5.0))
                self.b_input.setValue(lattice.get("b", 5.0))
                self.c_input.setValue(lattice.get("c", 5.0))
                self.alpha_input.setValue(lattice.get("alpha", 90.0))
                self.beta_input.setValue(lattice.get("beta", 90.0))
                self.gamma_input.setValue(lattice.get("gamma", 90.0))

            if "energy" in state:
                self.energy_input.setValue(state["energy"])

            if "file_path" in state and state["file_path"]:
                self.file_path_input.setText(state["file_path"])

            if "current_tab" in state:
                self.tab_widget.setCurrentIndex(state["current_tab"])

            return True
        except Exception:
            return False
