#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=no-name-in-module, import-error
from PyQt5.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QFormLayout,
    QGroupBox,
    QLabel,
    QPushButton,
    QDoubleSpinBox,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QHeaderView,
    QRadioButton,
    QButtonGroup,
    QFileDialog,
    QMessageBox,
)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QColor, QBrush
import csv
import os


class RangeInputWidget(QWidget):
    """Widget for input range (start, end, num_points)."""

    def __init__(self, label, parent=None):
        super().__init__(parent)

        # Main layout
        layout = QFormLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        # Create group with label
        group = QGroupBox(label)
        group_layout = QHBoxLayout(group)
        group_layout.setContentsMargins(10, 10, 10, 10)

        # Start value input with label
        start_widget = QWidget()
        start_layout = QFormLayout(start_widget)
        start_layout.setContentsMargins(0, 0, 0, 0)

        self.start_input = QDoubleSpinBox()
        self.start_input.setRange(-10.0, 10.0)
        self.start_input.setDecimals(3)
        self.start_input.setValue(0.0)
        start_layout.addRow("Start:", self.start_input)

        group_layout.addWidget(start_widget)

        # End value input with label
        end_widget = QWidget()
        end_layout = QFormLayout(end_widget)
        end_layout.setContentsMargins(0, 0, 0, 0)

        self.end_input = QDoubleSpinBox()
        self.end_input.setRange(-10.0, 10.0)
        self.end_input.setDecimals(3)
        self.end_input.setValue(-0.3)
        end_layout.addRow("End:", self.end_input)

        group_layout.addWidget(end_widget)

        layout.addRow(group)

    def get_range(self):
        """Get start and end values."""
        return (
            self.start_input.value(),
            self.end_input.value(),
        )

    def set_range(self, start, end):
        """Set start and end values."""
        self.start_input.setValue(start)
        self.end_input.setValue(end)

    def set_enabled(self, enabled):
        """Enable or disable widget."""
        self.start_input.setEnabled(enabled)
        self.end_input.setEnabled(enabled)


class HKLScanControls(QWidget):
    """Widget for HKL scan controls."""

    # Signal emitted when calculate button is clicked
    calculateClicked = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)

        # Main layout
        main_layout = QVBoxLayout(self)

        # Fixed tth input
        tth_group = QGroupBox("Fixed 2θ")
        tth_layout = QFormLayout(tth_group)

        self.tth_input = QDoubleSpinBox()
        self.tth_input.setRange(0.0, 180.0)
        self.tth_input.setValue(150.0)
        self.tth_input.setSuffix(" °")
        tth_layout.addRow("tth:", self.tth_input)

        main_layout.addWidget(tth_group)

        # Fixed angle selection (chi or phi)
        fixed_angle_group = QGroupBox("Fixed Sample Angle")
        fixed_angle_layout = QVBoxLayout(fixed_angle_group)

        angle_selection = QWidget()
        angle_selection_layout = QHBoxLayout(angle_selection)

        self.fix_chi_radio = QRadioButton("Fix χ")
        self.fix_phi_radio = QRadioButton("Fix φ")
        self.fix_chi_radio.setChecked(True)  # Default to fixed chi

        angle_selection_layout.addWidget(self.fix_chi_radio)
        angle_selection_layout.addWidget(self.fix_phi_radio)

        fixed_angle_layout.addWidget(angle_selection)

        # Create angle value inputs
        angle_values = QWidget()
        angle_values_layout = QHBoxLayout(angle_values)
        angle_values_layout.setContentsMargins(0, 0, 0, 0)

        # Chi input
        chi_widget = QWidget()
        chi_layout = QFormLayout(chi_widget)
        chi_layout.setContentsMargins(0, 0, 0, 0)
        self.chi_input = QDoubleSpinBox()
        self.chi_input.setRange(-180.0, 180.0)
        self.chi_input.setValue(0.0)
        self.chi_input.setSuffix(" °")
        chi_layout.addRow("χ:", self.chi_input)
        angle_values_layout.addWidget(chi_widget)

        # Phi input
        phi_widget = QWidget()
        phi_layout = QFormLayout(phi_widget)
        phi_layout.setContentsMargins(0, 0, 0, 0)
        self.phi_input = QDoubleSpinBox()
        self.phi_input.setRange(-180.0, 180.0)
        self.phi_input.setValue(0.0)
        self.phi_input.setSuffix(" °")
        phi_layout.addRow("φ:", self.phi_input)
        angle_values_layout.addWidget(phi_widget)

        fixed_angle_layout.addWidget(angle_values)
        main_layout.addWidget(fixed_angle_group)

        # HKL index selection
        hkl_group = QGroupBox("HKL Scan")
        hkl_layout = QVBoxLayout(hkl_group)

        # Index selection
        index_selection = QWidget()
        index_layout = QHBoxLayout(index_selection)
        index_layout.setContentsMargins(0, 0, 0, 0)

        self.h_toggle = QRadioButton("Deactivate H")
        self.k_toggle = QRadioButton("Deactivate K")
        self.l_toggle = QRadioButton("Deactivate L")
        self.l_toggle.setChecked(True)  # Default to deactivated L

        # Create a button group for mutual exclusion
        index_button_group = QButtonGroup(self)
        index_button_group.addButton(self.h_toggle)
        index_button_group.addButton(self.k_toggle)
        index_button_group.addButton(self.l_toggle)

        index_layout.addWidget(self.h_toggle)
        index_layout.addWidget(self.k_toggle)
        index_layout.addWidget(self.l_toggle)

        hkl_layout.addWidget(index_selection)

        # Create range widgets for H, K, L
        ranges_widget = QWidget()
        ranges_layout = QVBoxLayout(ranges_widget)
        ranges_layout.setContentsMargins(0, 0, 0, 0)

        self.h_range = RangeInputWidget("H Range")
        self.k_range = RangeInputWidget("K Range")
        self.l_range = RangeInputWidget("L Range")

        ranges_layout.addWidget(self.h_range)
        ranges_layout.addWidget(self.k_range)
        ranges_layout.addWidget(self.l_range)

        hkl_layout.addWidget(ranges_widget)

        # Number of points
        points_widget = QWidget()
        points_layout = QFormLayout(points_widget)
        points_layout.setContentsMargins(0, 0, 0, 0)

        self.num_points = QSpinBox()
        self.num_points.setRange(2, 100)
        self.num_points.setValue(10)
        points_layout.addRow("Number of points:", self.num_points)

        hkl_layout.addWidget(points_widget)

        main_layout.addWidget(hkl_group)

        # Calculate button
        self.calculate_button = QPushButton("Calculate")
        self.calculate_button.clicked.connect(self.calculateClicked.emit)
        main_layout.addWidget(self.calculate_button)

        # Connect signals
        self.h_toggle.toggled.connect(self._update_widget_states)
        self.k_toggle.toggled.connect(self._update_widget_states)
        self.l_toggle.toggled.connect(self._update_widget_states)
        self.fix_chi_radio.toggled.connect(self._update_fixed_angle_ui)
        self.fix_phi_radio.toggled.connect(self._update_fixed_angle_ui)

        # Initialize widget states
        self._update_widget_states()
        self._update_fixed_angle_ui()

    def _update_widget_states(self):
        """Update enabled state of widgets based on current selection."""
        self.h_range.set_enabled(not self.h_toggle.isChecked())
        self.k_range.set_enabled(not self.k_toggle.isChecked())
        self.l_range.set_enabled(not self.l_toggle.isChecked())

    def _update_fixed_angle_ui(self):
        """Update UI based on which angle is fixed."""
        is_chi_fixed = self.fix_chi_radio.isChecked()
        self.chi_input.setEnabled(is_chi_fixed)
        self.phi_input.setEnabled(not is_chi_fixed)

    def get_scan_parameters(self):
        """Get parameters for scan."""
        # Get deactivated index
        deactivated_index = None
        if self.h_toggle.isChecked():
            deactivated_index = "H"
        elif self.k_toggle.isChecked():
            deactivated_index = "K"
        else:  # l_toggle is checked
            deactivated_index = "L"

        # Get fixed angle
        fixed_angle_name = "chi" if self.fix_chi_radio.isChecked() else "phi"
        fixed_angle_value = (
            self.chi_input.value()
            if self.fix_chi_radio.isChecked()
            else self.phi_input.value()
        )

        # Get ranges
        h_start, h_end = self.h_range.get_range()
        k_start, k_end = self.k_range.get_range()
        l_start, l_end = self.l_range.get_range()

        return {
            "tth": self.tth_input.value(),
            "deactivated_index": deactivated_index,
            "fixed_angle_name": fixed_angle_name,
            "fixed_angle": fixed_angle_value,
            "start_points": (h_start, k_start, l_start),
            "end_points": (h_end, k_end, l_end),
            "num_points": self.num_points.value(),
        }


class HKLScanResultsTable(QTableWidget):
    """Table to display HKL scan results with multiple solutions."""

    def __init__(self, parent=None):
        super().__init__(parent)

        # Set up table
        self.setColumnCount(6)
        self.setHorizontalHeaderLabels(["H", "K", "L", "θ (°)", "φ (°)", "χ (°)"])
        self.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

        # Hide vertical header (row numbers)
        self.verticalHeader().setVisible(False)

        # Define colors for alternating groups white and gray
        self.group_colors = [
            QColor(255, 255, 255),  # White
            QColor(230, 230, 240),  # Gray
        ]

        # Enable sorting
        self.setSortingEnabled(True)

        # Add export button
        self.layout_wrapper = QVBoxLayout()
        self.layout_wrapper.setContentsMargins(0, 0, 0, 0)

        self.export_button = QPushButton("Export to CSV")
        self.export_button.clicked.connect(self.export_to_csv)
        self.export_button.setEnabled(False)  # Initially disabled until we have results

        self.layout_wrapper.addWidget(self)
        self.layout_wrapper.addWidget(self.export_button)

        # Store the last results for export
        self.last_results = None

    def display_results(self, results):
        """Display results in the table."""
        self.setSortingEnabled(False)  # Temporarily disable sorting
        self.setRowCount(0)  # Clear table

        # Check if we have results
        if not results or not results.get("success", False):
            self.export_button.setEnabled(False)
            self.last_results = None
            return

        # Store results for later export
        self.last_results = results

        # Get data from results
        h_values = results["H"]
        k_values = results["K"]
        l_values = results["L"]
        tth_values = results["tth"]  # Still get tth for export, just don't display it
        theta_values = results["theta"]
        phi_values = results["phi"]
        chi_values = results["chi"]

        # Group rows by common HK/KL/HL values
        deactivated_index = results.get("deactivated_index", None)
        groups = self._group_by_hkl(h_values, k_values, l_values, deactivated_index)

        # Color counter for alternating group colors
        current_color_idx = 0

        # Process each group and add to table with merged HKL values
        for group_indices in groups:
            if not group_indices:  # Skip empty groups
                continue

            # Select color for this group
            group_color = self.group_colors[current_color_idx % len(self.group_colors)]
            current_color_idx += 1

            # For each group, we should combine identical HKL values
            unique_hkls = {}  # Dictionary to track unique HKL combinations

            # First, identify unique HKL combinations within this group
            for i in group_indices:
                hkl_key = (
                    round(h_values[i], 6),
                    round(k_values[i], 6),
                    round(l_values[i], 6),
                )
                if hkl_key not in unique_hkls:
                    unique_hkls[hkl_key] = []
                unique_hkls[hkl_key].append(i)

            # Now add rows for each unique HKL combination
            for hkl_key, indices in unique_hkls.items():
                h, k, l = hkl_key
                first_row = self.rowCount()

                # Add a row for each angle solution for this HKL
                for j, i in enumerate(indices):
                    row_position = self.rowCount()
                    self.insertRow(row_position)

                    # For the first row of this HKL, add the HKL values
                    if j == 0:
                        h_item = QTableWidgetItem(f"{h:.4f}")
                        k_item = QTableWidgetItem(f"{k:.4f}")
                        l_item = QTableWidgetItem(f"{l:.4f}")

                        self.setItem(row_position, 0, h_item)
                        self.setItem(row_position, 1, k_item)
                        self.setItem(row_position, 2, l_item)

                    # Add angle values
                    self.setItem(
                        row_position, 3, QTableWidgetItem(f"{theta_values[i]:.1f}")
                    )
                    self.setItem(
                        row_position, 4, QTableWidgetItem(f"{phi_values[i]:.1f}")
                    )
                    self.setItem(
                        row_position, 5, QTableWidgetItem(f"{chi_values[i]:.1f}")
                    )

                    # Apply color to all cells in the row
                    for col in range(self.columnCount()):
                        item = self.item(row_position, col)
                        if item:
                            item.setBackground(QBrush(group_color))

                # If there are multiple solutions for this HKL, set row spans
                if len(indices) > 1:
                    self.setSpan(first_row, 0, len(indices), 1)  # H column
                    self.setSpan(first_row, 1, len(indices), 1)  # K column
                    self.setSpan(first_row, 2, len(indices), 1)  # L column

        # Re-enable sorting and export button
        self.setSortingEnabled(
            False
        )  # Keep sorting disabled as it interferes with row spans
        self.export_button.setEnabled(True)

    def _group_by_hkl(self, h_values, k_values, l_values, deactivated_index=None):
        """Group indices by common HK/KL/HL values."""
        groups = []

        # If deactivated_index is None, try to guess it based on constant values
        if deactivated_index is None:
            if len(set(h_values)) == 1:
                deactivated_index = "H"
            elif len(set(k_values)) == 1:
                deactivated_index = "K"
            elif len(set(l_values)) == 1:
                deactivated_index = "L"

        # If still None, default to grouping by all three
        if deactivated_index is None:
            # Group by all HKL values
            hkl_dict = {}
            for i in range(len(h_values)):
                key = (
                    round(h_values[i], 6),
                    round(k_values[i], 6),
                    round(l_values[i], 6),
                )
                if key not in hkl_dict:
                    hkl_dict[key] = []
                hkl_dict[key].append(i)

            # Convert dictionary to list of groups
            for indices in hkl_dict.values():
                groups.append(indices)
        else:
            # Group by the two indices that are not deactivated
            group_dict = {}
            for i in range(len(h_values)):
                if deactivated_index == "H":
                    key = (round(k_values[i], 6), round(l_values[i], 6))
                elif deactivated_index == "K":
                    key = (round(h_values[i], 6), round(l_values[i], 6))
                else:  # deactivated_index == 'L'
                    key = (round(h_values[i], 6), round(k_values[i], 6))

                if key not in group_dict:
                    group_dict[key] = []
                group_dict[key].append(i)

            # Convert dictionary to list of groups
            for indices in group_dict.values():
                groups.append(indices)

        return groups

    def export_to_csv(self):
        """Export results to a CSV file."""
        if not self.last_results:
            QMessageBox.warning(self, "Export Error", "No results to export.")
            return

        # Open file dialog to get save location
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save Results", "", "CSV Files (*.csv);;All Files (*)"
        )

        if not file_path:
            return  # User cancelled

        # Add .csv extension if not present
        if not file_path.endswith(".csv"):
            file_path += ".csv"

        try:
            with open(file_path, "w", newline="") as csvfile:
                writer = csv.writer(csvfile)

                # Write header
                writer.writerow(
                    [
                        "H",
                        "K",
                        "L",
                        "tth (deg)",
                        "theta (deg)",
                        "phi (deg)",
                        "chi (deg)",
                    ]
                )

                # Write data
                for i in range(len(self.last_results["tth"])):
                    writer.writerow(
                        [
                            f"{self.last_results['H'][i]:.6f}",
                            f"{self.last_results['K'][i]:.6f}",
                            f"{self.last_results['L'][i]:.6f}",
                            f"{self.last_results['tth'][i]:.6f}",
                            f"{self.last_results['theta'][i]:.6f}",
                            f"{self.last_results['phi'][i]:.6f}",
                            f"{self.last_results['chi'][i]:.6f}",
                        ]
                    )

            QMessageBox.information(
                self, "Export Success", f"Results exported to {file_path}"
            )
        except Exception as e:
            QMessageBox.critical(
                self, "Export Error", f"Error exporting results: {str(e)}"
            )

    def get_widget(self):
        """Return the widget containing the table and export button."""
        container = QWidget()
        container.setLayout(self.layout_wrapper)
        return container
