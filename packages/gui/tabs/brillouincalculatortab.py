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
    QTableWidget,
    QTableWidgetItem,
    QHeaderView,
    QFrame,
)
from PyQt5.QtCore import Qt, pyqtSlot, QMimeData
from PyQt5.QtGui import QDragEnterEvent, QDropEvent, QFont, QColor, QBrush
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
from packages.visualizer.unitcell_visualizer import UnitcellVisualizer
from packages.helpers.tips import Tips, set_tip
from packages.gui.components.hkl_scan_components import (
    HKLScanControls,
    HKLScanResultsTable,
)

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
        self.angles_to_hkl_visualizer = ScatteringVisualizer()  # first subtab
        self.hkl_to_angles_visualizer = ScatteringVisualizer()  # second subtab
        self.hk_fixed_tth_visualizer = ScatteringVisualizer()  # third subtab
        
        # Unit cell visualizers for each subtab
        self.angles_to_hkl_unitcell_viz = UnitcellVisualizer()
        self.hkl_to_angles_unitcell_viz = UnitcellVisualizer()
        self.hk_fixed_tth_unitcell_viz = UnitcellVisualizer()
        self.funtional_objects = [
            self.calculator,
            self.angles_to_hkl_visualizer,
            self.hkl_to_angles_visualizer,
            self.hk_fixed_tth_visualizer,
        ]  # group up the funtional objects and initailize them later on
        self.tips = Tips()

        # Store parameters for display
        self.parameters = None

        # Initialize UI
        super().__init__(main_window)
        self.init_ui()
        params = self.main_window.get_parameters()
        self.set_parameters(params)
        # Set window title
        self.setWindowTitle("Brillouin Zone Calculator")

    def init_ui(self):
        """Initialize UI components."""
        # Create tab widget for input methods
        self.tab_widget = QTabWidget()
        self.layout.addWidget(self.tab_widget, 0, 0)  # Back to (0,0) position

        # Create and add parameter header to the right
        self._create_parameter_header()

        # Create tabs for different functionalities
        self._create_angles_to_hkl_tab()
        self._create_hkl_to_angles_tab()
        self._create_hk_to_angles_tth_fixed_tab()
        self._create_hkl_scan_tab()  # Add the new HKL scan tab

    def _create_parameter_header(self):
        """Create parameter display header."""
        # Create main header frame
        header_frame = QFrame()
        header_frame.setObjectName("parameterPanel")
        header_frame.setFrameStyle(QFrame.StyledPanel | QFrame.Raised)
        header_frame.setLineWidth(1)
        header_frame.setFixedWidth(200)  # Reduced from 280 to 220 for narrower panel

        # Create header layout - vertical for sidebar
        header_layout = QVBoxLayout(header_frame)
        header_layout.setContentsMargins(15, 15, 15, 15)
        header_layout.setSpacing(15)  # Reduced from 25 to 15 for more compact layout

        # Crystal Structure section
        crystal_title = QLabel("Crystal Structure")
        crystal_title.setObjectName("parameterSectionTitle")
        header_layout.addWidget(crystal_title)

        self.crystal_info_label = QLabel(
            "a = -- Å\nb = -- Å\nc = -- Å\nα = --°\nβ = --°\nγ = --°"
        )
        self.crystal_info_label.setObjectName("parameterText")
        self.crystal_info_label.setWordWrap(True)
        header_layout.addWidget(self.crystal_info_label)

        # X-ray section
        xray_title = QLabel("X-ray Parameters")
        xray_title.setObjectName("parameterSectionTitle")
        header_layout.addWidget(xray_title)

        self.xray_info_label = QLabel("Energy: -- eV\nWavelength: -- Å")
        self.xray_info_label.setObjectName("parameterText")
        self.xray_info_label.setWordWrap(True)
        header_layout.addWidget(self.xray_info_label)

        # Add stretch to push edit button to the bottom
        header_layout.addStretch()

        # Edit button
        edit_button = QPushButton("Reset Parameters")
        edit_button.setObjectName("editParametersButton")
        edit_button.setToolTip("Return to parameter initialization window")
        edit_button.clicked.connect(self._edit_parameters)
        edit_button.setFixedHeight(35)  # Reduced from 45 to 35
        header_layout.addWidget(edit_button)

        # Add header to main layout - right side
        self.layout.addWidget(header_frame, 0, 1)

        # Set column stretch so tab widget takes most space
        self.layout.setColumnStretch(0, 5)  # Tab widget gets more space
        self.layout.setColumnStretch(1, 1)  # Header gets less space

    def set_parameters(self, params: dict):
        """Set parameters from global settings."""
        # Store parameters
        self.parameters = params

        # Update header display
        self._update_parameter_display()

        # Initialize functional objects
        for obj in self.funtional_objects:
            obj.initialize(params=params)

        # Initialize unit cell visualizers if CIF file is provided
        cif_file = params.get('cif_file')
        if cif_file:
            self._update_unitcell_visualizers(cif_file)

        print("params set!!!!!!!!")

    def _update_unitcell_visualizers(self, cif_file_path: str):
        """Update all unit cell visualizers with the CIF file."""
        try:
            unit_cell_vizs = [
                self.angles_to_hkl_unitcell_viz,
                self.hkl_to_angles_unitcell_viz,
                self.hk_fixed_tth_unitcell_viz
            ]
            
            for viz in unit_cell_vizs:
                viz.set_parameters({"cif_file": cif_file_path})
                viz.visualize_unitcell()
                
        except Exception as e:
            print(f"Error updating unit cell visualizers: {e}")

    def _update_parameter_display(self):
        """Update the parameter display in the header."""
        if not self.parameters:
            return

        # Update crystal structure display - vertical format
        crystal_text = (
            f"a = {self.parameters.get('a', 0):.2f} Å\n"
            f"b = {self.parameters.get('b', 0):.2f} Å\n"
            f"c = {self.parameters.get('c', 0):.2f} Å\n"
            f"α = {self.parameters.get('alpha', 0):.1f}°\n"
            f"β = {self.parameters.get('beta', 0):.1f}°\n"
            f"γ = {self.parameters.get('gamma', 0):.1f}°"
        )
        self.crystal_info_label.setText(crystal_text)

        # Update X-ray display - vertical format
        energy = self.parameters.get("energy", 0)
        # Convert eV to Angstrom: λ = hc/E = 12398.4 / E(eV)
        wavelength = 12398.4 / energy if energy > 0 else 0
        xray_text = f"Energy: {energy:.2f} eV\nWavelength: {wavelength:.3f} Å"
        self.xray_info_label.setText(xray_text)

    @pyqtSlot()
    def _edit_parameters(self):
        """Return to parameter initialization window."""
        if self.main_window:
            self.main_window.reset_parameters()

    def _set_tip(self, widget, name):
        """Set the tooltip and status tip for a widget by the name"""
        set_tip(widget, self.tips.tip(name))

    def _create_angles_to_hkl_tab(self):
        """Create tab for angles to HKL calculation."""
        angles_tab = QWidget()
        angles_layout = QHBoxLayout(angles_tab)

        # Left column - Parameters and inputs
        left_column = QWidget()
        left_layout = QVBoxLayout(left_column)

        # Input form
        form_group = QGroupBox("Scattering Angles")
        form_layout = QFormLayout(form_group)

        self.tth_angle_input = QDoubleSpinBox()
        self.tth_angle_input.setRange(0.0, 180.0)
        self.tth_angle_input.setValue(150.0)
        self.tth_angle_input.setSuffix(" °")
        self._set_tip(self.tth_angle_input, "TTH")
        form_layout.addRow("tth:", self.tth_angle_input)

        self.theta_angle_input = QDoubleSpinBox()
        self.theta_angle_input.setRange(-180.0, 180.0)
        self.theta_angle_input.setValue(50.0)
        self.theta_angle_input.setSuffix(" °")
        self._set_tip(self.theta_angle_input, "THETA")
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

        left_layout.addWidget(form_group)

        # Calculate button with fancy styling
        calculate_button = QPushButton("Calculate HKL")
        calculate_button.clicked.connect(self.calculate_hkl)
        calculate_button.setObjectName("calculateHKLButton")
        left_layout.addWidget(calculate_button)

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

        left_layout.addWidget(results_group)
        left_layout.addStretch()  # Add stretch to push content to top

        # Right column - Visualizers
        right_column = QWidget()
        right_layout = QVBoxLayout(right_column)

        # Scattering visualizer
        self.angles_to_hkl_visualizer.visualize_lab_system(is_clear=True)
        self.angles_to_hkl_visualizer.visualize_scattering_geometry(is_clear=False)
        right_layout.addWidget(self.angles_to_hkl_visualizer)
        
        # Unit cell visualizer
        right_layout.addWidget(self.angles_to_hkl_unitcell_viz)

        # Add columns to main layout
        angles_layout.addWidget(left_column, 1)  # Left column takes 1 part
        angles_layout.addWidget(right_column, 1.5)  # Right column takes 1.5 parts

        # Add to tab widget
        self.tab_widget.addTab(angles_tab, "Angles → HKL")

    def _create_hkl_to_angles_tab(self):
        """Create tab for HKL to angles calculation."""
        hkl_tab = QWidget()
        hkl_layout = QHBoxLayout(hkl_tab)

        # Left column - Parameters and inputs
        left_column = QWidget()
        left_layout = QVBoxLayout(left_column)

        # Input form
        form_group = QGroupBox("HKL Indices")
        form_layout = QFormLayout(form_group)

        self.H_input = QDoubleSpinBox()
        self.H_input.setRange(-10.0, 10.0)
        self.H_input.setDecimals(4)
        self.H_input.setValue(0.15)
        form_layout.addRow("H:", self.H_input)

        self.K_input = QDoubleSpinBox()
        self.K_input.setRange(-10.0, 10.0)
        self.K_input.setDecimals(4)
        self.K_input.setValue(0.1)
        form_layout.addRow("K:", self.K_input)

        self.L_input = QDoubleSpinBox()
        self.L_input.setRange(-10.0, 10.0)
        self.L_input.setDecimals(4)
        self.L_input.setValue(-0.5)
        form_layout.addRow("L:", self.L_input)

        left_layout.addWidget(form_group)

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
        self.fix_chi_radio.toggled.connect(self._update_fixed_angle_ui)
        self.fix_phi_radio.toggled.connect(self._update_fixed_angle_ui)

        # Initialize UI state
        self._update_fixed_angle_ui()

        left_layout.addWidget(constraints_group)

        # Calculate button with fancy styling
        calculate_button = QPushButton("Calculate Angles")
        calculate_button.clicked.connect(self.calculate_angles)
        calculate_button.setObjectName("calculateButton")
        left_layout.addWidget(calculate_button)

        # Results group
        results_group = QGroupBox("Results")
        results_layout = QVBoxLayout(results_group)

        # Add a table widget to display all possible solutions
        self.angles_results_table = QTableWidget()
        self.angles_results_table.setColumnCount(4)
        self.angles_results_table.setHorizontalHeaderLabels(
            ["tth (°)", "θ (°)", "φ (°)", "χ (°)"]
        )
        self.angles_results_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.Stretch
        )
        # Connect selection change signal
        self.angles_results_table.itemSelectionChanged.connect(
            self.on_angle_solution_selected
        )
        results_layout.addWidget(self.angles_results_table)

        # Add clear button with smaller styling
        clear_button = QPushButton("Clear Results")
        clear_button.clicked.connect(self.clear_hkl_to_angles_results)
        clear_button.setObjectName("clearButton")
        results_layout.addWidget(clear_button)

        left_layout.addWidget(results_group)
        left_layout.addStretch()  # Add stretch to push content to top

        # Right column - Visualizers
        right_column = QWidget()
        right_layout = QVBoxLayout(right_column)

        # Scattering visualizer
        self.hkl_to_angles_visualizer.visualize_lab_system(is_clear=True)
        self.hkl_to_angles_visualizer.visualize_scattering_geometry(is_clear=False)
        right_layout.addWidget(self.hkl_to_angles_visualizer)
        
        # Unit cell visualizer
        right_layout.addWidget(self.hkl_to_angles_unitcell_viz)

        # Add columns to main layout
        hkl_layout.addWidget(left_column, 1)  # Left column takes 1 part
        hkl_layout.addWidget(right_column, 1.5)  # Right column takes 1.5 parts

        # Add to tab widget
        self.tab_widget.addTab(hkl_tab, "HKL → Angles")

    def _create_hk_to_angles_tth_fixed_tab(self):
        """Create tab for HK to angles calculation with fixed tth."""
        hk_tab = QWidget()
        hk_layout = QHBoxLayout(hk_tab)

        # Left column - Parameters and inputs
        left_column = QWidget()
        left_layout = QVBoxLayout(left_column)

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
        self.chi_toggle.toggled.connect(self._update_fixed_angle_ui_tth)
        self.phi_toggle.toggled.connect(self._update_fixed_angle_ui_tth)

        # Ensure only one angle is deactivated at a time
        self.chi_toggle.toggled.connect(
            lambda checked: self._ensure_one_angle_deactivated("chi", checked)
        )
        self.phi_toggle.toggled.connect(
            lambda checked: self._ensure_one_angle_deactivated("phi", checked)
        )

        # Initialize UI state
        self._update_fixed_angle_ui_tth()

        left_layout.addWidget(constraints_group)

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
        self.H_input_tth.setDecimals(4)
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
        self.K_input_tth.setDecimals(4)
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
        self.L_input_tth.setDecimals(4)
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
        self.H_toggle.toggled.connect(self._update_free_index_ui)
        self.K_toggle.toggled.connect(self._update_free_index_ui)
        self.L_toggle.toggled.connect(self._update_free_index_ui)

        # Ensure only one index is free at a time
        self.H_toggle.toggled.connect(
            lambda checked: self._ensure_one_free_index("H", checked)
        )
        self.K_toggle.toggled.connect(
            lambda checked: self._ensure_one_free_index("K", checked)
        )
        self.L_toggle.toggled.connect(
            lambda checked: self._ensure_one_free_index("L", checked)
        )

        left_layout.addWidget(hkl_group)

        # Calculate button with fancy styling
        calculate_button = QPushButton("Calculate Angles")
        calculate_button.clicked.connect(self.calculate_angles_tth_fixed)
        calculate_button.setObjectName("calculateButton")
        left_layout.addWidget(calculate_button)

        # Results group
        results_group = QGroupBox("Results")
        results_layout = QVBoxLayout(results_group)

        # Add a table widget to display all possible solutions
        self.angles_results_table_tth = QTableWidget()
        self.angles_results_table_tth.setColumnCount(4)
        self.angles_results_table_tth.setHorizontalHeaderLabels(
            ["tth (°)", "θ (°)", "φ (°)", "χ (°)"]
        )
        self.angles_results_table_tth.horizontalHeader().setSectionResizeMode(
            QHeaderView.Stretch
        )
        # Connect selection change signal
        self.angles_results_table_tth.itemSelectionChanged.connect(
            self.on_angle_solution_selected_tth
        )
        results_layout.addWidget(self.angles_results_table_tth)

        # Add clear button with smaller styling
        clear_button_tth = QPushButton("Clear Results")
        clear_button_tth.clicked.connect(self.clear_hk_tth_fixed_results)
        clear_button_tth.setObjectName("clearButton")
        results_layout.addWidget(clear_button_tth)

        left_layout.addWidget(results_group)
        left_layout.addStretch()  # Add stretch to push content to top

        # Right column - Visualizers
        right_column = QWidget()
        right_layout = QVBoxLayout(right_column)

        # Scattering visualizer
        self.hk_fixed_tth_visualizer.visualize_lab_system(is_clear=True)
        self.hk_fixed_tth_visualizer.visualize_scattering_geometry(is_clear=False)
        right_layout.addWidget(self.hk_fixed_tth_visualizer)
        
        # Unit cell visualizer
        right_layout.addWidget(self.hk_fixed_tth_unitcell_viz)

        # Add columns to main layout
        hk_layout.addWidget(left_column, 1)  # Left column takes 1 part
        hk_layout.addWidget(right_column, 1.5)  # Right column takes 1.5 parts

        # Add to tab widget
        self.tab_widget.addTab(hk_tab, "HK to Angles | tth fixed")

    def _create_hkl_scan_tab(self):
        """Create tab for scanning a range of HKL values."""
        scan_tab = QWidget()
        scan_layout = QHBoxLayout(scan_tab)

        # Create controls widget
        self.hkl_scan_controls = HKLScanControls()
        self.hkl_scan_controls.calculateClicked.connect(self.calculate_hkl_scan)
        scan_layout.addWidget(self.hkl_scan_controls, 1)

        # Create results table
        self.hkl_scan_results_table = HKLScanResultsTable()
        scan_layout.addWidget(
            self.hkl_scan_results_table.get_widget(), 3.3
        )  # Increased from 3 to 3.3 for 10% wider

        # Add to tab widget
        self.tab_widget.addTab(scan_tab, "HKL Scan | tth fixed")

    @pyqtSlot()
    def calculate_hkl_scan(self):
        """Calculate angles for a range of HKL values."""
        try:
            # Check if calculator is initialized
            if not self.calculator.is_initialized():
                QMessageBox.warning(
                    self, "Warning", "Please initialize the calculator first!"
                )
                self.tab_widget.setCurrentIndex(0)
                return

            # Get parameters for scan
            params = self.hkl_scan_controls.get_scan_parameters()

            # Calculate angles for the scan
            result = self.calculator.calculate_angles_tth_fixed_scan(
                tth=params["tth"],
                start_points=params["start_points"],
                end_points=params["end_points"],
                num_points=params["num_points"],
                deactivated_index=params["deactivated_index"],
                fixed_angle_name=params["fixed_angle_name"],
                fixed_angle=params["fixed_angle"],
            )

            # Check for success
            success = result.get("success", False)
            if not success:
                QMessageBox.warning(
                    self, "Warning", result.get("error", "Unknown error")
                )
                return

            # Display results in table
            self.hkl_scan_results_table.display_results(result)

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error calculating HKL scan: {str(e)}")

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
            self.angles_to_hkl_visualizer.visualize_lab_system(
                is_clear=True, chi=chi, phi=phi
            )
            self.angles_to_hkl_visualizer.visualize_scattering_geometry(
                scattering_angles=result, is_clear=False
            )

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error calculating HKL: {str(e)}")

    @pyqtSlot()
    def _update_fixed_angle_ui(self):
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
                H=self.H_input.value(),
                K=self.K_input.value(),
                L=self.L_input.value(),
                fixed_angle=fixed_angle_value,
                fixed_angle_name=fixed_angle_name,
            )

            success = result.get("success", False)
            if not success:
                QMessageBox.warning(
                    self, "Warning", result.get("error", "Unknown error")
                )
                return

            # Extract results - these are now lists
            tth_values = result["tth"]
            theta_values = result["theta"]
            phi_values = result["phi"]
            chi_values = result["chi"]

            # Store the current row count to know where the new results start
            start_row = self.angles_results_table.rowCount()
            num_solutions = len(tth_values)

            # First, clear highlighting from all existing rows (set to white background)
            for row in range(start_row):
                for col in range(self.angles_results_table.columnCount()):
                    item = self.angles_results_table.item(row, col)
                    if item:
                        item.setBackground(QBrush(QColor(255, 255, 255)))  # White background

            # Append new solutions to the table
            for i in range(num_solutions):
                row_position = self.angles_results_table.rowCount()
                self.angles_results_table.insertRow(row_position)

                # Add items to the row with highlighting
                tth_item = QTableWidgetItem(f"{tth_values[i]:.1f}")
                theta_item = QTableWidgetItem(f"{theta_values[i]:.1f}")
                phi_item = QTableWidgetItem(f"{phi_values[i]:.1f}")
                chi_item = QTableWidgetItem(f"{chi_values[i]:.1f}")
                
                # Set background color for new items (230, 230, 240)
                highlight_color = QBrush(QColor(230, 230, 240))
                tth_item.setBackground(highlight_color)
                theta_item.setBackground(highlight_color)
                phi_item.setBackground(highlight_color)
                chi_item.setBackground(highlight_color)

                self.angles_results_table.setItem(row_position, 0, tth_item)
                self.angles_results_table.setItem(row_position, 1, theta_item)
                self.angles_results_table.setItem(row_position, 2, phi_item)
                self.angles_results_table.setItem(row_position, 3, chi_item)

            # Scroll to the bottom to show the new results
            self.angles_results_table.scrollToBottom()

            # If we have at least one solution, visualize the first one
            if num_solutions > 0:
                first_solution = {
                    "tth": tth_values[0],
                    "theta": theta_values[0],
                    "phi": phi_values[0],
                    "chi": chi_values[0],
                    "H": result["H"],
                    "K": result["K"],
                    "L": result["L"],
                }

                # Update visualization with the first solution
                self.hkl_to_angles_visualizer.visualize_lab_system(
                    is_clear=True, chi=chi_values[0], phi=phi_values[0]
                )
                self.hkl_to_angles_visualizer.visualize_scattering_geometry(
                    scattering_angles=first_solution, is_clear=False
                )

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error calculating angles: {str(e)}")

    @pyqtSlot()
    def _update_fixed_angle_ui_tth(self):
        """Update UI based on which angle is fixed."""
        self.chi_input_tth.setEnabled(not self.chi_toggle.isChecked())
        self.phi_input_tth.setEnabled(not self.phi_toggle.isChecked())

    @pyqtSlot()
    def _update_free_index_ui(self):
        """Update UI based on which index is free."""
        self.H_input_tth.setEnabled(not self.H_toggle.isChecked())
        self.K_input_tth.setEnabled(not self.K_toggle.isChecked())
        self.L_input_tth.setEnabled(not self.L_toggle.isChecked())

    @pyqtSlot()
    def _ensure_one_free_index(self, index, checked):
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
            fixed_angle_name = "phi" if self.chi_toggle.isChecked() else "chi"
            fixed_angle_value = (
                self.phi_input_tth.value()
                if self.chi_toggle.isChecked()
                else self.chi_input_tth.value()
            )

            # Get input values
            tth = self.tth_fixed_input.value()
            H = self.H_input_tth.value() if fixed_index != "H" else None
            K = self.K_input_tth.value() if fixed_index != "K" else None
            L = self.L_input_tth.value() if fixed_index != "L" else None

            # Calculate angles
            result = self.calculator.calculate_angles_tth_fixed(
                tth=tth,
                H=H,
                K=K,
                L=L,
                fixed_angle_name=fixed_angle_name,
                fixed_angle=fixed_angle_value,
            )

            success = result.get("success", False)
            if not success:
                QMessageBox.warning(
                    self, "Warning", result.get("error", "Unknown error")
                )
                return

            # Extract results - these are now lists
            tth_values = result["tth"]
            theta_values = result["theta"]
            phi_values = result["phi"]
            chi_values = result["chi"]

            # Store the current row count to know where the new results start
            start_row = self.angles_results_table_tth.rowCount()
            num_solutions = len(tth_values)

            # First, clear highlighting from all existing rows (set to white background)
            for row in range(start_row):
                for col in range(self.angles_results_table_tth.columnCount()):
                    item = self.angles_results_table_tth.item(row, col)
                    if item:
                        item.setBackground(QBrush(QColor(255, 255, 255)))  # White background

            # Append new solutions to the table
            for i in range(num_solutions):
                row_position = self.angles_results_table_tth.rowCount()
                self.angles_results_table_tth.insertRow(row_position)

                # Add items to the row with highlighting
                tth_item = QTableWidgetItem(f"{tth_values[i]:.1f}")
                theta_item = QTableWidgetItem(f"{theta_values[i]:.1f}")
                phi_item = QTableWidgetItem(f"{phi_values[i]:.1f}")
                chi_item = QTableWidgetItem(f"{chi_values[i]:.1f}")
                
                # Set background color for new items (230, 230, 240)
                highlight_color = QBrush(QColor(230, 230, 240))
                tth_item.setBackground(highlight_color)
                theta_item.setBackground(highlight_color)
                phi_item.setBackground(highlight_color)
                chi_item.setBackground(highlight_color)

                self.angles_results_table_tth.setItem(row_position, 0, tth_item)
                self.angles_results_table_tth.setItem(row_position, 1, theta_item)
                self.angles_results_table_tth.setItem(row_position, 2, phi_item)
                self.angles_results_table_tth.setItem(row_position, 3, chi_item)

            # Scroll to the bottom to show the new results
            self.angles_results_table_tth.scrollToBottom()

            # If we have at least one solution, update the HKL values and visualize it
            if num_solutions > 0:
                # Update the HKL values
                H_result = result["H"]
                K_result = result["K"]
                L_result = result["L"]

                # Update the input fields without triggering signals
                self.H_input_tth.blockSignals(True)
                self.K_input_tth.blockSignals(True)
                self.L_input_tth.blockSignals(True)

                self.H_input_tth.setValue(H_result)
                self.K_input_tth.setValue(K_result)
                self.L_input_tth.setValue(L_result)

                self.H_input_tth.blockSignals(False)
                self.K_input_tth.blockSignals(False)
                self.L_input_tth.blockSignals(False)

                # Create a first solution for visualization
                first_solution = {
                    "tth": tth_values[0],
                    "theta": theta_values[0],
                    "phi": phi_values[0],
                    "chi": chi_values[0],
                    "H": H_result,
                    "K": K_result,
                    "L": L_result,
                }

                # Update visualization with the first solution
                self.hk_fixed_tth_visualizer.visualize_lab_system(
                    is_clear=True, chi=chi_values[0], phi=phi_values[0]
                )
                self.hk_fixed_tth_visualizer.visualize_scattering_geometry(
                    scattering_angles=first_solution, is_clear=False
                )

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error calculating angles: {str(e)}")

    @pyqtSlot()
    def _ensure_one_angle_deactivated(self, angle, checked):
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

    @pyqtSlot()
    def on_angle_solution_selected(self):
        """Handle selection of a specific angle solution from the results table."""
        selected_items = self.angles_results_table.selectedItems()
        if not selected_items:
            return

        # Get the row of the first selected item
        row = selected_items[0].row()

        # Get values from the selected row - convert from text to float
        tth = float(self.angles_results_table.item(row, 0).text())
        theta = float(self.angles_results_table.item(row, 1).text())
        phi = float(self.angles_results_table.item(row, 2).text())
        chi = float(self.angles_results_table.item(row, 3).text())

        # Create a solution dictionary
        solution = {
            "tth": tth,
            "theta": theta,
            "phi": phi,
            "chi": chi,
            "H": self.H_input.value(),
            "K": self.K_input.value(),
            "L": self.L_input.value(),
        }

        # Update visualization with the selected solution
        self.hkl_to_angles_visualizer.visualize_lab_system(
            is_clear=True, chi=chi, phi=phi
        )
        self.hkl_to_angles_visualizer.visualize_scattering_geometry(
            scattering_angles=solution, is_clear=False
        )

    @pyqtSlot()
    def on_angle_solution_selected_tth(self):
        """Handle selection of a specific angle solution from the results table."""
        selected_items = self.angles_results_table_tth.selectedItems()
        if not selected_items:
            return

        # Get the row of the first selected item
        row = selected_items[0].row()

        # Get values from the selected row - convert from text to float
        tth = float(self.angles_results_table_tth.item(row, 0).text())
        theta = float(self.angles_results_table_tth.item(row, 1).text())
        phi = float(self.angles_results_table_tth.item(row, 2).text())
        chi = float(self.angles_results_table_tth.item(row, 3).text())

        # Create a solution dictionary
        solution = {
            "tth": tth,
            "theta": theta,
            "phi": phi,
            "chi": chi,
            "H": self.H_input_tth.value(),
            "K": self.K_input_tth.value(),
            "L": self.L_input_tth.value(),
        }

        # Update visualization with the selected solution
        self.hk_fixed_tth_visualizer.visualize_lab_system(
            is_clear=True, chi=chi, phi=phi
        )
        self.hk_fixed_tth_visualizer.visualize_scattering_geometry(
            scattering_angles=solution, is_clear=False
        )

    @pyqtSlot()
    def clear_hkl_to_angles_results(self):
        """Clear all results from the HKL to angles table."""
        self.angles_results_table.setRowCount(0)

    @pyqtSlot()
    def clear_hk_tth_fixed_results(self):
        """Clear all results from the HK to angles (tth fixed) table."""
        self.angles_results_table_tth.setRowCount(0)

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
