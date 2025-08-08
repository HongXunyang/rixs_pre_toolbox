#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=no-name-in-module, import-error
from PyQt5.QtWidgets import (
    QGridLayout,
    QFormLayout,
    QLabel,
    QPushButton,
    QDoubleSpinBox,
    QGroupBox,
    QMessageBox,
    QVBoxLayout,
    QTableWidget,
    QTableWidgetItem,
    QHeaderView,
    QTextEdit,
)
from PyQt5.QtCore import pyqtSlot
import sys
import os
import numpy as np

# Add parent directory to path to allow imports from sibling packages
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from packages.gui.tabs.tab_interface import TabInterface
from packages.structure_factor_calculator.interface import StructureFactorCalculator
from packages.visualizer.structure_factor_visualizer import StructureFactorVisualizer
from packages.helpers.tips import Tips, set_tip


class StructureFactorTab(TabInterface):
    """Tab for calculating structure factors using X-ray scattering."""

    def __init__(self, main_window=None):
        # Create backend instance
        self.calculator = StructureFactorCalculator()
        self.visualizer = StructureFactorVisualizer()
        self.tips = Tips()

        # Track initialization state
        self._calculator_initialized = False
        # CIF path is provided globally via InitWindow parameters

        # Initialize UI first
        super().__init__(main_window)
        self.setWindowTitle("Structure Factor Calculator")

    def init_ui(self):
        """Initialize UI components."""
        self._create_structure_factors_ui()

    def set_parameters(self, params: dict):
        """Set parameters from global lattice configuration.

        Note: Structure factor calculator requires CIF file and energy,
        which are not part of the global lattice parameters.
        """
        # The visualizer can be initialized with basic parameters
        self.visualizer.initialize(params)

        # For the calculator, we need CIF file and energy to be set separately
        # This will be handled through the UI controls

    def _set_tip(self, widget, name):
        """Set the tooltip and status tip for a widget by the name"""
        set_tip(widget, self.tips.tip(name))

    def _create_structure_factors_ui(self):
        """Create UI for calculating structure factors."""
        main_layout = QGridLayout()
        self.layout.addLayout(main_layout, 0, 0)

        # Configuration section (energy + initialize)
        config_group = QGroupBox("Configuration")
        config_layout = QFormLayout(config_group)

        # Energy input (in eV)
        self.energy_input = QDoubleSpinBox()
        self.energy_input.setRange(1.0, 100000.0)  # 1 eV to 100 keV
        self.energy_input.setDecimals(1)
        self.energy_input.setValue(10000.0)  # Default 10 keV
        self.energy_input.setSuffix(" eV")
        config_layout.addRow("X-ray Energy:", self.energy_input)

        # Initialize button
        self.init_btn = QPushButton("Initialize Calculator")
        self.init_btn.clicked.connect(self.initialize_calculator)
        config_layout.addRow("", self.init_btn)

        # Status display
        self.status_label = QLabel(
            "Status: Provide CIF in initialization window, then initialize"
        )
        self.status_label.setStyleSheet("color: orange; font-weight: bold;")
        config_layout.addRow("", self.status_label)

        main_layout.addWidget(config_group, 0, 0)

        # HKL input section
        hkl_group = QGroupBox("HKL Input")
        hkl_layout = QVBoxLayout(hkl_group)

        # Option 1: Use default HKL list
        self.use_default_btn = QPushButton("Use Default HKL List")
        self.use_default_btn.clicked.connect(self.use_default_hkl)
        self.use_default_btn.setEnabled(False)
        hkl_layout.addWidget(self.use_default_btn)

        # Option 2: Custom HKL input
        custom_hkl_layout = QFormLayout()
        self.custom_hkl_input = QTextEdit()
        self.custom_hkl_input.setMaximumHeight(100)
        self.custom_hkl_input.setPlaceholderText(
            "Enter HKL values, one per line: h k l\nExample:\n1 0 0\n0 1 0\n1 1 1"
        )
        custom_hkl_layout.addRow("Custom HKL:", self.custom_hkl_input)

        self.use_custom_btn = QPushButton("Use Custom HKL List")
        self.use_custom_btn.clicked.connect(self.use_custom_hkl)
        self.use_custom_btn.setEnabled(False)
        custom_hkl_layout.addRow("", self.use_custom_btn)

        hkl_layout.addLayout(custom_hkl_layout)
        main_layout.addWidget(hkl_group, 1, 0)

        # Results section
        results_group = QGroupBox("Results")
        results_layout = QVBoxLayout(results_group)

        self.results_table = QTableWidget()
        self.results_table.setColumnCount(4)
        self.results_table.setHorizontalHeaderLabels(
            ["H", "K", "L", "Structure Factor |F|"]
        )
        self.results_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        results_layout.addWidget(self.results_table)

        main_layout.addWidget(results_group, 2, 0)

        # Visualization panel on the right
        self.plot_group = QGroupBox("Visualization")
        plot_layout = QVBoxLayout(self.plot_group)
        plot_layout.addWidget(self.visualizer)
        main_layout.addWidget(self.plot_group, 0, 1, 3, 1)

    @pyqtSlot()
    def initialize_calculator(self):
        """Initialize the structure factor calculator."""
        try:
            # Get CIF path from global parameters provided by InitWindow
            params = self.main_window.get_parameters() if self.main_window else None
            cif_path = params.get("cif_file") if params else None
            if not cif_path:
                QMessageBox.warning(
                    self,
                    "Missing CIF",
                    "Please load a valid CIF file in the initialization window first.",
                )
                return

            energy = self.energy_input.value()

            # Initialize the calculator with CIF file and energy
            self.calculator.initialize(cif_path, energy)

            self._calculator_initialized = True
            self.status_label.setText("Status: Calculator initialized successfully")
            self.status_label.setStyleSheet("color: green; font-weight: bold;")

            # Enable HKL input buttons
            self.use_default_btn.setEnabled(True)
            self.use_custom_btn.setEnabled(True)

        except Exception as e:
            QMessageBox.critical(
                self, "Error", f"Failed to initialize calculator: {str(e)}"
            )
            self.status_label.setText(f"Status: Initialization failed - {str(e)}")
            self.status_label.setStyleSheet("color: red; font-weight: bold;")

    @pyqtSlot()
    def use_default_hkl(self):
        """Calculate structure factors using default HKL list."""
        if not self._calculator_initialized:
            QMessageBox.warning(
                self, "Warning", "Please initialize the calculator first!"
            )
            return

        try:
            # Use the default HKL list from the interface
            results = self.calculator.calculate_structure_factors()

            # Get the default HKL list
            from packages.structure_factor_calculator.interface import DEFAULT_HKL_LIST

            hkl_list = DEFAULT_HKL_LIST

            # Debug output
            print(f"Default HKL list: {hkl_list}")
            print(f"Results shape: {np.array(results).shape}")
            print(f"Results sample: {results[:3] if len(results) > 0 else 'empty'}")
            print(
                f"Absolute values: {np.abs(results)[:3] if len(results) > 0 else 'empty'}"
            )

            # Update results table and visualization
            self.update_results_table(hkl_list, results)
            success = self.visualizer.visualize_structure_factors(
                hkl_list, np.abs(results)
            )
            print(f"Visualization success: {success}")

        except Exception as e:
            QMessageBox.critical(
                self, "Error", f"Error calculating structure factors: {str(e)}"
            )

    @pyqtSlot()
    def use_custom_hkl(self):
        """Calculate structure factors using custom HKL list."""
        if not self._calculator_initialized:
            QMessageBox.warning(
                self, "Warning", "Please initialize the calculator first!"
            )
            return

        try:
            # Parse custom HKL input
            hkl_list = self._parse_custom_hkl()
            if not hkl_list:
                QMessageBox.warning(self, "Warning", "Please enter valid HKL values!")
                return

            # Calculate structure factors
            results = self.calculator.calculate_structure_factors(hkl_list)

            # Update results table and visualization
            self.update_results_table(hkl_list, results)
            self.visualizer.visualize_structure_factors(hkl_list, np.abs(results))

        except Exception as e:
            QMessageBox.critical(
                self, "Error", f"Error calculating structure factors: {str(e)}"
            )

    def _parse_custom_hkl(self):
        """Parse custom HKL input from text widget."""
        text = self.custom_hkl_input.toPlainText().strip()
        if not text:
            return []

        hkl_list = []
        try:
            for line in text.split("\n"):
                line = line.strip()
                if line and not line.startswith("#"):  # Skip empty lines and comments
                    parts = line.split()
                    if len(parts) >= 3:
                        h, k, l = int(parts[0]), int(parts[1]), int(parts[2])
                        hkl_list.append([h, k, l])
        except ValueError as e:
            QMessageBox.warning(self, "Warning", f"Invalid HKL format: {str(e)}")
            return []

        return hkl_list

    def update_results_table(self, hkl_list, sf_values):
        """Update the results table with calculated structure factors."""
        # Clear the table
        self.results_table.setRowCount(0)

        # Add results to the table
        for i, (hkl, sf) in enumerate(zip(hkl_list, sf_values)):
            h, k, l = hkl[0], hkl[1], hkl[2]
            sf_magnitude = np.abs(sf)  # Get magnitude of complex structure factor

            self.results_table.insertRow(i)
            self.results_table.setItem(i, 0, QTableWidgetItem(str(h)))
            self.results_table.setItem(i, 1, QTableWidgetItem(str(k)))
            self.results_table.setItem(i, 2, QTableWidgetItem(str(l)))
            self.results_table.setItem(i, 3, QTableWidgetItem(f"{sf_magnitude:.4f}"))

    def get_module_instance(self):
        """Get the backend module instance."""
        return self.calculator

    def clear(self):
        """Clear all inputs and results."""
        self.custom_hkl_input.clear()
        self.results_table.setRowCount(0)
        self.energy_input.setValue(10000.0)
        self._calculator_initialized = False
        self.status_label.setText(
            "Status: Provide CIF in initialization window, then initialize"
        )
        self.status_label.setStyleSheet("color: orange; font-weight: bold;")
        self.use_default_btn.setEnabled(False)
        self.use_custom_btn.setEnabled(False)

        # Clear visualization
        self.visualizer.clear_plot()

    def get_state(self):
        """Get the current state for session saving."""
        return {
            "energy": self.energy_input.value(),
            "custom_hkl": self.custom_hkl_input.toPlainText(),
            "calculator_initialized": self._calculator_initialized,
        }

    def set_state(self, state):
        """Restore tab state from saved session."""
        try:
            if "energy" in state:
                self.energy_input.setValue(state["energy"])

            if "custom_hkl" in state:
                self.custom_hkl_input.setPlainText(state["custom_hkl"])

            if "calculator_initialized" in state and state["calculator_initialized"]:
                # Try to reinitialize if we have the required data globally
                params = self.main_window.get_parameters() if self.main_window else None
                if (
                    params
                    and params.get("cif_file")
                    and os.path.exists(params.get("cif_file"))
                ):
                    self.initialize_calculator()

            return True
        except Exception:
            return False
