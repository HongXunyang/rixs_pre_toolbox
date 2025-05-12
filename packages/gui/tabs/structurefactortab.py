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
    QFileDialog,
    QMessageBox,
    QVBoxLayout,
    QHBoxLayout,
    QTableWidget,
    QTableWidgetItem,
    QHeaderView,
)
from PyQt5.QtCore import Qt, pyqtSlot
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
    """Tab for calculating structure factors."""

    def __init__(self, main_window=None):
        # Create backend instance
        self.calculator = StructureFactorCalculator()
        self.visualizer = StructureFactorVisualizer()
        self.funtional_objects = [self.calculator, self.visualizer]
        self.tips = Tips()

        # Initialize UI
        super().__init__(main_window)
        self.setWindowTitle("Structure Factor Calculator")

        # Set parameters if the main window has them
        params = self.main_window.get_parameters() if self.main_window else None
        if params:
            self.set_parameters(params)

    def init_ui(self):
        """Initialize UI components."""
        # Create the structure factors UI directly in the main layout
        self._create_structure_factors_ui()

    def set_parameters(self, params: dict):
        """Set parameters from global settings."""
        for obj in self.funtional_objects:
            obj.initialize(params=params)

    def _set_tip(self, widget, name):
        """Set the tooltip and status tip for a widget by the name"""
        set_tip(widget, self.tips.tip(name))

    def _create_structure_factors_ui(self):
        """Create UI for calculating structure factors."""
        main_layout = QGridLayout()
        self.layout.addLayout(main_layout, 0, 0)

        # Input form
        input_group = QGroupBox("HKL Range")
        input_layout = QFormLayout(input_group)

        # H range
        h_range_layout = QHBoxLayout()
        self.h_min_input = QDoubleSpinBox()
        self.h_min_input.setRange(-10.0, 10.0)
        self.h_min_input.setDecimals(0)
        self.h_min_input.setValue(-3)
        h_range_layout.addWidget(self.h_min_input)

        h_range_layout.addWidget(QLabel("to"))

        self.h_max_input = QDoubleSpinBox()
        self.h_max_input.setRange(-10.0, 10.0)
        self.h_max_input.setDecimals(0)
        self.h_max_input.setValue(3)
        h_range_layout.addWidget(self.h_max_input)
        input_layout.addRow("H range:", h_range_layout)

        # K range
        k_range_layout = QHBoxLayout()
        self.k_min_input = QDoubleSpinBox()
        self.k_min_input.setRange(-10.0, 10.0)
        self.k_min_input.setDecimals(0)
        self.k_min_input.setValue(-3)
        k_range_layout.addWidget(self.k_min_input)

        k_range_layout.addWidget(QLabel("to"))

        self.k_max_input = QDoubleSpinBox()
        self.k_max_input.setRange(-10.0, 10.0)
        self.k_max_input.setDecimals(0)
        self.k_max_input.setValue(3)
        k_range_layout.addWidget(self.k_max_input)
        input_layout.addRow("K range:", k_range_layout)

        # L range
        l_range_layout = QHBoxLayout()
        self.l_min_input = QDoubleSpinBox()
        self.l_min_input.setRange(-10.0, 10.0)
        self.l_min_input.setDecimals(0)
        self.l_min_input.setValue(-3)
        l_range_layout.addWidget(self.l_min_input)

        l_range_layout.addWidget(QLabel("to"))

        self.l_max_input = QDoubleSpinBox()
        self.l_max_input.setRange(-10.0, 10.0)
        self.l_max_input.setDecimals(0)
        self.l_max_input.setValue(3)
        l_range_layout.addWidget(self.l_max_input)
        input_layout.addRow("L range:", l_range_layout)

        main_layout.addWidget(input_group, 0, 0)

        # Calculate button
        calculate_button = QPushButton("Calculate Structure Factors")
        calculate_button.clicked.connect(self.calculate_structure_factors)
        main_layout.addWidget(calculate_button, 1, 0)

        # Results table
        results_group = QGroupBox("Results")
        results_layout = QVBoxLayout(results_group)

        self.results_table = QTableWidget()
        self.results_table.setColumnCount(5)  # Added Intensity column
        self.results_table.setHorizontalHeaderLabels(
            ["H", "K", "L", "Structure Factor", "Intensity"]
        )
        self.results_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        results_layout.addWidget(self.results_table)

        main_layout.addWidget(results_group, 2, 0)

        # Visualization panel on the right
        self.plot_group = QGroupBox("Visualization")
        plot_layout = QVBoxLayout(self.plot_group)

        # Use the structure factor visualizer
        plot_layout.addWidget(self.visualizer)

        main_layout.addWidget(self.plot_group, 0, 1, 3, 1)

    @pyqtSlot()
    def calculate_structure_factors(self):
        """Calculate structure factors for multiple HKL values."""
        try:
            # Check if calculator is initialized
            if not self.calculator.is_initialized():
                QMessageBox.warning(
                    self, "Warning", "Please initialize the calculator first!"
                )
                return

            # Get HKL ranges
            h_range = range(
                int(self.h_min_input.value()), int(self.h_max_input.value()) + 1
            )
            k_range = range(
                int(self.k_min_input.value()), int(self.k_max_input.value()) + 1
            )
            l_range = range(
                int(self.l_min_input.value()), int(self.l_max_input.value()) + 1
            )

            # Prepare HKL list
            hkl_list = []
            for h in h_range:
                for k in k_range:
                    for l in l_range:
                        hkl_list.append((h, k, l))

            # Calculate structure factors using the backend
            results = self.calculator.calculate_structure_factors(hkl_list)

            # If the backend doesn't implement the calculation yet, use a placeholder
            if results[0] is None:
                # Placeholder implementation until backend is fully implemented
                results = [
                    abs(h + k + l) for h, k, l in hkl_list
                ]  # Simple placeholder formula

            # Update results table
            self.update_results_table(hkl_list, results)

            # Update visualization using the structure factor visualizer
            self.visualizer.visualize_structure_factors(hkl_list, results)

        except Exception as e:
            QMessageBox.critical(
                self, "Error", f"Error calculating structure factors: {str(e)}"
            )

    def update_results_table(self, hkl_list, sf_values):
        """Update the results table with calculated structure factors."""
        # Clear the table
        self.results_table.setRowCount(0)

        # Add results to the table
        for i, ((h, k, l), sf) in enumerate(zip(hkl_list, sf_values)):
            intensity = sf**2  # Calculate intensity as |F|²
            self.results_table.insertRow(i)
            self.results_table.setItem(i, 0, QTableWidgetItem(str(h)))
            self.results_table.setItem(i, 1, QTableWidgetItem(str(k)))
            self.results_table.setItem(i, 2, QTableWidgetItem(str(l)))
            self.results_table.setItem(i, 3, QTableWidgetItem(f"{sf:.4f}"))
            self.results_table.setItem(i, 4, QTableWidgetItem(f"{intensity:.4f}"))

    def get_module_instance(self):
        """Get the backend module instance."""
        return self.calculator

    def get_state(self):
        """Get the current state for session saving."""
        return {
            "hkl_range": {
                "h_min": self.h_min_input.value(),
                "h_max": self.h_max_input.value(),
                "k_min": self.k_min_input.value(),
                "k_max": self.k_max_input.value(),
                "l_min": self.l_min_input.value(),
                "l_max": self.l_max_input.value(),
            }
        }

    def set_state(self, state):
        """Restore tab state from saved session."""
        try:
            if "hkl_range" in state:
                hkl_range = state["hkl_range"]
                self.h_min_input.setValue(hkl_range.get("h_min", -3))
                self.h_max_input.setValue(hkl_range.get("h_max", 3))
                self.k_min_input.setValue(hkl_range.get("k_min", -3))
                self.k_max_input.setValue(hkl_range.get("k_max", 3))
                self.l_min_input.setValue(hkl_range.get("l_min", -3))
                self.l_max_input.setValue(hkl_range.get("l_max", 3))

            return True
        except Exception:
            return False
