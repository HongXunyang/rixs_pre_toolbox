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
        self.create_hk_to_angles_tth_fixed_tab()

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

        # Fix angle selection
        fix_angle_group = QGroupBox("Fix Angle")
        fix_angle_layout = QHBoxLayout(fix_angle_group)

        self.fix_chi_radio = QRadioButton("Fix χ")
        self.fix_phi_radio = QRadioButton("Fix φ")
        self.fix_chi_radio.setChecked(True)  # Default to fixed chi

        fix_angle_layout.addWidget(self.fix_chi_radio)
        fix_angle_layout.addWidget(self.fix_phi_radio)

        constraints_layout.addRow(fix_angle_group)

        # Create a horizontal layout for chi and phi inputs
        angles_row = QWidget()
        angles_layout = QHBoxLayout(angles_row)
        angles_layout.setContentsMargins(0, 0, 0, 0)

        # Chi input
        chi_widget = QWidget()
        chi_layout = QFormLayout(chi_widget)
        chi_layout.setContentsMargins(0, 0, 0, 0)
        self.chi_input = QDoubleSpinBox()
        self.chi_input.setRange(-180.0, 180.0)
        self.chi_input.setValue(0.0)
        self.chi_input.setSuffix(" °")
        chi_layout.addRow("χ:", self.chi_input)
        angles_layout.addWidget(chi_widget)

        # Phi input
        phi_widget = QWidget()
        phi_layout = QFormLayout(phi_widget)
        phi_layout.setContentsMargins(0, 0, 0, 0)
        self.phi_input = QDoubleSpinBox()
        self.phi_input.setRange(-180.0, 180.0)
        self.phi_input.setValue(0.0)
        self.phi_input.setSuffix(" °")
        phi_layout.addRow("φ:", self.phi_input)
        angles_layout.addWidget(phi_widget)

        constraints_layout.addRow(angles_row)

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

    def create_hk_to_angles_tth_fixed_tab(self):
        """Create tab for HK to angles calculation with fixed tth."""
        hk_tab = QWidget()
        hk_layout = QGridLayout(hk_tab)

        # Constraints group
        constraints_group = QGroupBox("Constraints")
        constraints_layout = QFormLayout(constraints_group)

        # tth input
        self.tth_fixed_input = QDoubleSpinBox()
        self.tth_fixed_input.setRange(0.0, 180.0)
        self.tth_fixed_input.setValue(150.0)
        self.tth_fixed_input.setSuffix(" °")
        constraints_layout.addRow("tth:", self.tth_fixed_input)

        # Create a vertical layout for chi and phi inputs
        angles_widget = QWidget()
        angles_layout = QVBoxLayout(angles_widget)
        angles_layout.setContentsMargins(0, 0, 0, 0)

        # Chi input row
        chi_row = QWidget()
        chi_layout = QHBoxLayout(chi_row)
        chi_layout.setContentsMargins(0, 0, 0, 0)
        chi_form = QWidget()
        chi_form_layout = QFormLayout(chi_form)
        chi_form_layout.setContentsMargins(0, 0, 0, 0)
        self.chi_input_tth = QDoubleSpinBox()
        self.chi_input_tth.setRange(-180.0, 180.0)
        self.chi_input_tth.setValue(0.0)
        self.chi_input_tth.setSuffix(" °")
        self.chi_input_tth.setEnabled(False)  # Default to disabled
        chi_form_layout.addRow("χ:", self.chi_input_tth)
        chi_layout.addWidget(chi_form)
        self.chi_toggle = QPushButton("De-activate")
        self.chi_toggle.setCheckable(True)
        self.chi_toggle.setChecked(True)  # Default to deactivated
        chi_layout.addWidget(self.chi_toggle)
        angles_layout.addWidget(chi_row)

        # Phi input row
        phi_row = QWidget()
        phi_layout = QHBoxLayout(phi_row)
        phi_layout.setContentsMargins(0, 0, 0, 0)
        phi_form = QWidget()
        phi_form_layout = QFormLayout(phi_form)
        phi_form_layout.setContentsMargins(0, 0, 0, 0)
        self.phi_input_tth = QDoubleSpinBox()
        self.phi_input_tth.setRange(-180.0, 180.0)
        self.phi_input_tth.setValue(0.0)
        self.phi_input_tth.setSuffix(" °")
        self.phi_input_tth.setEnabled(True)  # Default to enabled
        phi_form_layout.addRow("φ:", self.phi_input_tth)
        phi_layout.addWidget(phi_form)
        self.phi_toggle = QPushButton("De-activate")
        self.phi_toggle.setCheckable(True)
        self.phi_toggle.setChecked(False)  # Default to activated
        phi_layout.addWidget(self.phi_toggle)
        angles_layout.addWidget(phi_row)

        constraints_layout.addRow(angles_widget)

        # Connect toggle buttons to enable/disable corresponding inputs
        self.chi_toggle.toggled.connect(self.update_fixed_angle_ui_tth)
        self.phi_toggle.toggled.connect(self.update_fixed_angle_ui_tth)

        # Ensure only one angle is deactivated at a time
        self.chi_toggle.toggled.connect(
            lambda checked: self.ensure_one_angle_deactivated("chi", checked)
        )
        self.phi_toggle.toggled.connect(
            lambda checked: self.ensure_one_angle_deactivated("phi", checked)
        )

        # Initialize UI state
        self.update_fixed_angle_ui_tth()

        hk_layout.addWidget(constraints_group, 0, 0)

        # HKL indices group
        hkl_group = QGroupBox("HKL Indices")
        hkl_layout = QFormLayout(hkl_group)

        # Create a vertical layout for HKL inputs
        hkl_row = QWidget()
        hkl_input_layout = QVBoxLayout(hkl_row)
        hkl_input_layout.setContentsMargins(0, 0, 0, 0)

        # H input row
        h_row = QWidget()
        h_layout = QHBoxLayout(h_row)
        h_layout.setContentsMargins(0, 0, 0, 0)
        h_form = QWidget()
        h_form_layout = QFormLayout(h_form)
        h_form_layout.setContentsMargins(0, 0, 0, 0)
        self.H_input_tth = QDoubleSpinBox()
        self.H_input_tth.setRange(-10.0, 10.0)
        self.H_input_tth.setDecimals(3)
        self.H_input_tth.setValue(0.15)
        h_form_layout.addRow("H:", self.H_input_tth)
        h_layout.addWidget(h_form)
        self.H_toggle = QPushButton("De-activate")
        self.H_toggle.setCheckable(True)
        self.H_toggle.setChecked(False)  # Default to active (gray)
        h_layout.addWidget(self.H_toggle)
        hkl_input_layout.addWidget(h_row)

        # K input row
        k_row = QWidget()
        k_layout = QHBoxLayout(k_row)
        k_layout.setContentsMargins(0, 0, 0, 0)
        k_form = QWidget()
        k_form_layout = QFormLayout(k_form)
        k_form_layout.setContentsMargins(0, 0, 0, 0)
        self.K_input_tth = QDoubleSpinBox()
        self.K_input_tth.setRange(-10.0, 10.0)
        self.K_input_tth.setDecimals(3)
        self.K_input_tth.setValue(0.1)
        k_form_layout.addRow("K:", self.K_input_tth)
        k_layout.addWidget(k_form)
        self.K_toggle = QPushButton("De-activate")
        self.K_toggle.setCheckable(True)
        self.K_toggle.setChecked(False)  # Default to active (gray)
        k_layout.addWidget(self.K_toggle)
        hkl_input_layout.addWidget(k_row)

        # L input row
        l_row = QWidget()
        l_layout = QHBoxLayout(l_row)
        l_layout.setContentsMargins(0, 0, 0, 0)
        l_form = QWidget()
        l_form_layout = QFormLayout(l_form)
        l_form_layout.setContentsMargins(0, 0, 0, 0)
        self.L_input_tth = QDoubleSpinBox()
        self.L_input_tth.setRange(-10.0, 10.0)
        self.L_input_tth.setDecimals(3)
        self.L_input_tth.setValue(-0.5)
        self.L_input_tth.setEnabled(False)  # Default to disabled (grayed out)
        l_form_layout.addRow("L:", self.L_input_tth)
        l_layout.addWidget(l_form)
        self.L_toggle = QPushButton("De-activate")
        self.L_toggle.setCheckable(True)
        self.L_toggle.setChecked(True)  # Default to deactivated (red)
        l_layout.addWidget(self.L_toggle)
        hkl_input_layout.addWidget(l_row)

        hkl_layout.addRow(hkl_row)

        # Connect toggle buttons to enable/disable corresponding inputs
        self.H_toggle.toggled.connect(self.update_free_index_ui)
        self.K_toggle.toggled.connect(self.update_free_index_ui)
        self.L_toggle.toggled.connect(self.update_free_index_ui)

        # Ensure only one index is free at a time
        self.H_toggle.toggled.connect(
            lambda checked: self.ensure_one_free_index("H", checked)
        )
        self.K_toggle.toggled.connect(
            lambda checked: self.ensure_one_free_index("K", checked)
        )
        self.L_toggle.toggled.connect(
            lambda checked: self.ensure_one_free_index("L", checked)
        )

        hk_layout.addWidget(hkl_group, 1, 0)

        # Calculate button
        calculate_button = QPushButton("Calculate Angles")
        calculate_button.clicked.connect(self.calculate_angles_tth_fixed)
        hk_layout.addWidget(calculate_button, 2, 0)

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
        self.tth_result_tth = QLineEdit()
        self.set_tip(self.tth_result_tth, "TTH")
        self.tth_result_tth.setReadOnly(True)
        tth_layout.addRow("tth:", self.tth_result_tth)
        first_row_layout.addWidget(tth_widget)

        phi_widget = QWidget()
        phi_layout = QFormLayout(phi_widget)
        phi_layout.setContentsMargins(0, 0, 0, 0)
        self.phi_result_tth = QLineEdit()
        self.phi_result_tth.setReadOnly(True)
        phi_layout.addRow("φ:", self.phi_result_tth)
        first_row_layout.addWidget(phi_widget)

        results_layout.addWidget(first_row)

        # Second row: theta and chi
        second_row = QWidget()
        second_row_layout = QHBoxLayout(second_row)
        second_row_layout.setContentsMargins(0, 0, 0, 0)

        theta_widget = QWidget()
        theta_layout = QFormLayout(theta_widget)
        theta_layout.setContentsMargins(0, 0, 0, 0)
        self.theta_result_tth = QLineEdit()
        self.set_tip(self.theta_result_tth, "THETA")
        self.theta_result_tth.setReadOnly(True)
        theta_layout.addRow("θ:", self.theta_result_tth)
        second_row_layout.addWidget(theta_widget)

        chi_widget = QWidget()
        chi_layout = QFormLayout(chi_widget)
        chi_layout.setContentsMargins(0, 0, 0, 0)
        self.chi_result_tth = QLineEdit()
        self.chi_result_tth.setReadOnly(True)
        chi_layout.addRow("χ:", self.chi_result_tth)
        second_row_layout.addWidget(chi_widget)

        results_layout.addWidget(second_row)

        hk_layout.addWidget(results_group, 3, 0)

        # Visualizer
        self.hk_visualizer = ScatteringVisualizer()
        self.hk_visualizer.visualize_lab_system(is_clear=True)
        self.hk_visualizer.visualize_scattering_geometry(is_clear=False)
        hk_layout.addWidget(self.hk_visualizer, 0, 1, 3, 1)

        # Add to tab widget
        self.tab_widget.addTab(hk_tab, "HK to Angles | tth fixed")

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

    @pyqtSlot()
    def update_fixed_angle_ui_tth(self):
        """Update UI based on which angle is fixed."""
        self.chi_input_tth.setEnabled(not self.chi_toggle.isChecked())
        self.phi_input_tth.setEnabled(not self.phi_toggle.isChecked())

    @pyqtSlot()
    def update_free_index_ui(self):
        """Update UI based on which index is free."""
        self.H_input_tth.setEnabled(not self.H_toggle.isChecked())
        self.K_input_tth.setEnabled(not self.K_toggle.isChecked())
        self.L_input_tth.setEnabled(not self.L_toggle.isChecked())

    @pyqtSlot()
    def ensure_one_free_index(self, index, checked):
        """Ensure only one index is free at a time."""
        if checked:
            # If this index is now free, uncheck the others
            if index != "H":
                self.H_toggle.setChecked(False)
            if index != "K":
                self.K_toggle.setChecked(False)
            if index != "L":
                self.L_toggle.setChecked(False)
        else:
            # If this index is now fixed, ensure at least one other is free
            if not (
                self.H_toggle.isChecked()
                or self.K_toggle.isChecked()
                or self.L_toggle.isChecked()
            ):
                # If no index is free, re-check this one
                if index == "H":
                    self.H_toggle.setChecked(True)
                elif index == "K":
                    self.K_toggle.setChecked(True)
                else:
                    self.L_toggle.setChecked(True)

    @pyqtSlot()
    def calculate_angles_tth_fixed(self):
        """Calculate angles from HK with fixed tth."""
        try:
            # Check if calculator is initialized
            if not self.calculator.is_initialized():
                QMessageBox.warning(
                    self, "Warning", "Please initialize the calculator first!"
                )
                self.tab_widget.setCurrentIndex(0)
                return

            # Determine which index to fix
            fixed_index = None
            if self.H_toggle.isChecked():
                fixed_index = "H"
            elif self.K_toggle.isChecked():
                fixed_index = "K"
            else:  # L_toggle is checked
                fixed_index = "L"

            # Determine which angle to fix
            fixed_angle_name = "chi" if self.chi_toggle.isChecked() else "phi"
            fixed_angle_value = (
                self.chi_input_tth.value()
                if self.chi_toggle.isChecked()
                else self.phi_input_tth.value()
            )

            # Get input values
            tth = self.tth_fixed_input.value()
            H = self.H_input_tth.value() if fixed_index != "H" else None
            K = self.K_input_tth.value() if fixed_index != "K" else None
            L = self.L_input_tth.value() if fixed_index != "L" else None

            # Calculate angles
            result = self.calculator.calculate_angles_tth_fixed(
                tth=tth,
                H_crystal=H,
                K_crystal=K,
                L_crystal=L,
                fixed_angle_name=fixed_angle_name,
                fixed_angle=fixed_angle_value,
            )

            success = result.get("success", False)
            if not success:
                QMessageBox.warning(
                    self, "Warning", result.get("error", "Unknown error")
                )
                return

            # Update results
            self.tth_result_tth.setText(f"{result['tth']:.4f}°")
            self.theta_result_tth.setText(f"{result['theta']:.4f}°")
            self.phi_result_tth.setText(f"{result['phi']:.4f}°")
            self.chi_result_tth.setText(f"{result['chi']:.4f}°")
            # update the H, K, L values
            self.H_input_tth.setValue(result["H"])
            self.K_input_tth.setValue(result["K"])
            self.L_input_tth.setValue(result["L"])
            # Update visualization
            self.hk_visualizer.visualize_lab_system(
                is_clear=True, chi=result["chi"], phi=result["phi"]
            )
            self.hk_visualizer.visualize_scattering_geometry(
                scattering_angles=result, is_clear=False
            )

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error calculating angles: {str(e)}")

    @pyqtSlot()
    def ensure_one_angle_deactivated(self, angle, checked):
        """Ensure only one angle is deactivated at a time."""
        if checked:
            # If this angle is now deactivated, activate the other
            if angle == "chi":
                self.phi_toggle.setChecked(False)
            else:
                self.chi_toggle.setChecked(False)
        else:
            # If this angle is now activated, ensure at least one other is deactivated
            if not (self.chi_toggle.isChecked() or self.phi_toggle.isChecked()):
                # If no angle is deactivated, deactivate this one
                if angle == "chi":
                    self.chi_toggle.setChecked(True)
                else:
                    self.phi_toggle.setChecked(True)

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
