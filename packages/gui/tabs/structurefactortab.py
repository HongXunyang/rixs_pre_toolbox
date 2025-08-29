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
    QTabWidget,
    QWidget,
    QSlider,
    QHBoxLayout,
    QSpinBox,
    QLineEdit,
    QStackedLayout,
)
from PyQt5.QtCore import Qt, pyqtSlot
import sys
import os
import numpy as np

# Add parent directory to path to allow imports from sibling packages
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from .tab_interface import TabInterface
from packages.structure_factor_calculator import StructureFactorCalculator
from packages.visualizer import StructureFactorVisualizer3D, StructureFactorVisualizer2D
from packages.helpers.tips import Tips, set_tip


class StructureFactorTab(TabInterface):
    """Tab for calculating structure factors using X-ray scattering."""

    def __init__(self, main_window=None):
        # Create backend instance
        self.calculator = StructureFactorCalculator()
        self.tips = Tips()

        # Visualizers for new In-use subtab
        self.inuse_visualizer3d = StructureFactorVisualizer3D()
        self.hk_plane_visualizer = StructureFactorVisualizer2D()
        self.hl_plane_visualizer = StructureFactorVisualizer2D()
        self.kl_plane_visualizer = StructureFactorVisualizer2D()

        # Visualizer for deprecated subtab
        self.visualizer3d = StructureFactorVisualizer3D()

        # Track initialization state
        self._calculator_initialized = False
        # CIF path is provided globally via InitWindow parameters

        # Initialize UI first
        super().__init__(main_window)
        self.setWindowTitle("Structure Factor Calculator")

    def init_ui(self):
        """Initialize UI components with subtabs."""
        self.tab_widget = QTabWidget()
        self.layout.addWidget(self.tab_widget, 0, 0)

        self._create_in_use_tab()
        self._create_customized_tab()

    def set_parameters(self, params: dict):
        """Set parameters from global lattice configuration.

        Note: Structure factor calculator requires CIF file and energy,
        which are not part of the global lattice parameters.
        """
        # Initialize visualizers
        self.inuse_visualizer3d.initialize(params)
        self.hk_plane_visualizer.clear_plot()
        self.hl_plane_visualizer.clear_plot()
        self.kl_plane_visualizer.clear_plot()

        # For the calculator, we need CIF file and energy to be set separately
        # This will be handled through the UI controls

    def _set_tip(self, widget, name):
        """Set the tooltip and status tip for a widget by the name"""
        set_tip(widget, self.tips.tip(name))

    def _create_in_use_tab(self):
        """Create the new primary subtab with 2x4 layout and live-updating 2D slices."""
        inuse_tab = QWidget()
        main_layout = QGridLayout(inuse_tab)

        # Row 0 left: 3D visualizer
        main_layout.addWidget(self.inuse_visualizer3d, 0, 0)

        # Plane toggle buttons in configuration section for compactness
        self.hk_toggle_btn = QPushButton("HK plane")
        self.hl_toggle_btn = QPushButton("HL plane")
        self.kl_toggle_btn = QPushButton("KL plane")
        for btn in (self.hk_toggle_btn, self.hl_toggle_btn, self.kl_toggle_btn):
            btn.setCheckable(True)
        self.hk_toggle_btn.clicked.connect(lambda: self._set_active_inuse_plane("HK"))
        self.hl_toggle_btn.clicked.connect(lambda: self._set_active_inuse_plane("HL"))
        self.kl_toggle_btn.clicked.connect(lambda: self._set_active_inuse_plane("KL"))

        # (1,0) Configuration section (energy + initialize)
        config_group = QGroupBox("Configuration")
        config_layout = QFormLayout(config_group)

        self.energy_input_inuse = QDoubleSpinBox()
        self.energy_input_inuse.setRange(1.0, 100000.0)
        self.energy_input_inuse.setDecimals(1)
        self.energy_input_inuse.setValue(10000.0)
        self.energy_input_inuse.setSuffix(" eV")
        config_layout.addRow("X-ray Energy:", self.energy_input_inuse)

        # Add plane toggle row
        toggle_row = QWidget()
        toggle_layout = QHBoxLayout(toggle_row)
        toggle_layout.setContentsMargins(0, 0, 0, 0)
        toggle_layout.addWidget(self.hk_toggle_btn)
        toggle_layout.addWidget(self.hl_toggle_btn)
        toggle_layout.addWidget(self.kl_toggle_btn)
        config_layout.addRow("Plane:", toggle_row)

        self.init_btn_inuse = QPushButton("Initialize Calculator")
        self.init_btn_inuse.clicked.connect(self.initialize_calculator_inuse)
        config_layout.addRow("", self.init_btn_inuse)

        self.status_label_inuse = QLabel(
            "Status: Provide CIF in initialization window, then initialize"
        )
        self.status_label_inuse.setStyleSheet("color: orange; font-weight: bold;")
        config_layout.addRow("", self.status_label_inuse)
        # Row 1 left: configuration group
        main_layout.addWidget(config_group, 1, 0)

        # Convenience ranges for integer grid
        self._grid_max = 5  # inclusive max for varying integer indices

        # Controls for HK plane (L fixed)
        self.hk_ctrl_group = self._create_fixed_index_controls(
            label_prefix="HK plane", fixed_name="L", default_value=0
        )
        self.hk_ctrl_group["spin"].valueChanged.connect(
            lambda v: self._sync_slider_and_update(
                self.hk_ctrl_group, v, self.update_hk_plane_plot
            )
        )
        # Update 3D overlay plane when control changes
        self.hk_ctrl_group["spin"].valueChanged.connect(
            lambda v: self.inuse_visualizer3d.set_plane_values(L=int(v))
        )
        self.hk_ctrl_group["slider"].valueChanged.connect(
            lambda v: self.inuse_visualizer3d.set_plane_values(L=int(v))
        )
        self.hk_ctrl_group["slider"].valueChanged.connect(
            lambda v: self._sync_spin_and_update(
                self.hk_ctrl_group, v, self.update_hk_plane_plot
            )
        )

        # Controls for HL plane (K fixed)
        self.hl_ctrl_group = self._create_fixed_index_controls(
            label_prefix="HL plane", fixed_name="K", default_value=0
        )
        self.hl_ctrl_group["spin"].valueChanged.connect(
            lambda v: self._sync_slider_and_update(
                self.hl_ctrl_group, v, self.update_hl_plane_plot
            )
        )
        self.hl_ctrl_group["spin"].valueChanged.connect(
            lambda v: self.inuse_visualizer3d.set_plane_values(K=int(v))
        )
        self.hl_ctrl_group["slider"].valueChanged.connect(
            lambda v: self._sync_spin_and_update(
                self.hl_ctrl_group, v, self.update_hl_plane_plot
            )
        )
        self.hl_ctrl_group["slider"].valueChanged.connect(
            lambda v: self.inuse_visualizer3d.set_plane_values(K=int(v))
        )

        # Controls for KL plane (H fixed)
        self.kl_ctrl_group = self._create_fixed_index_controls(
            label_prefix="KL plane", fixed_name="H", default_value=0
        )
        self.kl_ctrl_group["spin"].valueChanged.connect(
            lambda v: self._sync_slider_and_update(
                self.kl_ctrl_group, v, self.update_kl_plane_plot
            )
        )
        self.kl_ctrl_group["spin"].valueChanged.connect(
            lambda v: self.inuse_visualizer3d.set_plane_values(H=int(v))
        )
        self.kl_ctrl_group["slider"].valueChanged.connect(
            lambda v: self._sync_spin_and_update(
                self.kl_ctrl_group, v, self.update_kl_plane_plot
            )
        )
        self.kl_ctrl_group["slider"].valueChanged.connect(
            lambda v: self.inuse_visualizer3d.set_plane_values(H=int(v))
        )

        # --- Right column: use stacked layouts for 2D plane and its controls ---
        # Plane stack (top-right)
        self._plane_stack_container = QWidget()
        self._plane_stack = QStackedLayout(self._plane_stack_container)
        self._plane_stack.addWidget(self.hk_plane_visualizer)
        self._plane_stack.addWidget(self.hl_plane_visualizer)
        self._plane_stack.addWidget(self.kl_plane_visualizer)
        main_layout.addWidget(self._plane_stack_container, 0, 1)

        # Control stack (bottom-right)
        self._ctrl_stack_container = QWidget()
        self._ctrl_stack = QStackedLayout(self._ctrl_stack_container)
        self._ctrl_stack.addWidget(self.hk_ctrl_group["group"])  # index 0
        self._ctrl_stack.addWidget(self.hl_ctrl_group["group"])  # index 1
        self._ctrl_stack.addWidget(self.kl_ctrl_group["group"])  # index 2
        main_layout.addWidget(self._ctrl_stack_container, 1, 1)

        # Make columns 50/50 and give more vertical space to plots
        main_layout.setColumnStretch(0, 1)
        main_layout.setColumnStretch(1, 1)
        main_layout.setRowStretch(0, 3)
        main_layout.setRowStretch(1, 1)

        # Add the subtab as the first tab
        self.tab_widget.addTab(inuse_tab, "HKL plane")

        # Default active plane: HK
        self._set_active_inuse_plane("HK")

    def _update_toggle_styles(self, active: str):
        """Update toggle button colors based on active plane."""
        active_css = "background-color: #2ecc71; color: white; font-weight: bold;"
        inactive_css = "background-color: #bdc3c7; color: #333333;"
        mapping = {
            "HK": self.hk_toggle_btn,
            "HL": self.hl_toggle_btn,
            "KL": self.kl_toggle_btn,
        }
        for name, btn in mapping.items():
            if name == active:
                btn.setChecked(True)
                btn.setStyleSheet(active_css)
            else:
                btn.setChecked(False)
                btn.setStyleSheet(inactive_css)

    def _set_active_inuse_plane(self, plane: str):
        """Show only the selected 2D plane and its control panel."""
        plane = plane.upper()
        self._update_toggle_styles(plane)
        # Switch plane/control stacks instead of show/hide
        index_map = {"HK": 0, "HL": 1, "KL": 2}
        idx = index_map.get(plane, 0)
        self._plane_stack.setCurrentIndex(idx)
        self._ctrl_stack.setCurrentIndex(idx)
        # Update 3D plane highlighting
        if plane == "HK":
            self.inuse_visualizer3d.set_active_plane("L")  # HK plane means L constant
        elif plane == "HL":
            self.inuse_visualizer3d.set_active_plane("K")
        elif plane == "KL":
            self.inuse_visualizer3d.set_active_plane("H")
        # Update corresponding plot if initialized
        if self.calculator.is_initialized:
            if plane == "HK":
                self.update_hk_plane_plot()
            elif plane == "HL":
                self.update_hl_plane_plot()
            elif plane == "KL":
                self.update_kl_plane_plot()

    def _create_fixed_index_controls(
        self, label_prefix: str, fixed_name: str, default_value: int
    ):
        """Utility to create a group with an integer spinbox and slider tied together."""
        group = QGroupBox(f"{label_prefix}")
        layout = QFormLayout(group)
        row = QWidget()
        row_layout = QHBoxLayout(row)
        row_layout.setContentsMargins(0, 0, 0, 0)

        spin = QSpinBox()
        spin.setRange(0, 5)
        spin.setValue(default_value)
        slider = QSlider()
        slider.setOrientation(Qt.Horizontal)
        slider.setRange(0, 5)
        slider.setValue(default_value)

        row_layout.addWidget(spin)
        row_layout.addWidget(slider)
        layout.addRow(f"{fixed_name}:", row)

        return {"group": group, "spin": spin, "slider": slider}

    def _sync_slider_and_update(self, ctrl_group, value: int, updater):
        ctrl_group["slider"].blockSignals(True)
        ctrl_group["slider"].setValue(value)
        ctrl_group["slider"].blockSignals(False)
        updater()

    def _sync_spin_and_update(self, ctrl_group, value: int, updater):
        ctrl_group["spin"].blockSignals(True)
        ctrl_group["spin"].setValue(value)
        ctrl_group["spin"].blockSignals(False)
        updater()

    def _generate_plane_points(
        self, varying_a: str, varying_b: str, fixed_name: str, fixed_value: int
    ):
        """Create integer HKL points for a plane with two varying indices in [0, _grid_max]."""
        points = []
        for a in range(0, self._grid_max + 1):
            for b in range(0, self._grid_max + 1):
                values = {"H": 0, "K": 0, "L": 0}
                values[varying_a] = a
                values[varying_b] = b
                values[fixed_name] = fixed_value
                points.append([values["H"], values["K"], values["L"]])
        return points

    def _generate_hkl_cube(self, max_index: int = 5):
        """Generate a full integer HKL grid from 0..max_index for 3D visualization."""
        cube = []
        for h in range(0, max_index + 1):
            for k in range(0, max_index + 1):
                for l in range(0, max_index + 1):
                    cube.append([h, k, l])
        return cube

    def update_hk_plane_plot(self):
        if not self.calculator.is_initialized:
            return
        L_val = int(self.hk_ctrl_group["spin"].value())
        hkl_list = self._generate_plane_points("H", "K", "L", L_val)
        results = self.calculator.calculate_structure_factors(hkl_list)
        # compute |F(0,0,L_val)| as color scale max reference
        ref = self.calculator.calculate_structure_factors([[0, 0, 0]])
        value_max = float(np.abs(ref[0])) if len(ref) > 0 else None
        arr = np.array(hkl_list)
        self.hk_plane_visualizer.visualize_plane(
            arr[:, 0], arr[:, 1], np.abs(results), "H", "K", "L", L_val, value_max
        )

    def update_hl_plane_plot(self):
        if not self.calculator.is_initialized:
            return
        K_val = int(self.hl_ctrl_group["spin"].value())
        hkl_list = self._generate_plane_points("H", "L", "K", K_val)
        results = self.calculator.calculate_structure_factors(hkl_list)
        # compute |F(0,K_val,0)| as color scale max reference
        ref = self.calculator.calculate_structure_factors([[0, 0, 0]])
        value_max = float(np.abs(ref[0])) if len(ref) > 0 else None
        arr = np.array(hkl_list)
        self.hl_plane_visualizer.visualize_plane(
            arr[:, 0], arr[:, 2], np.abs(results), "H", "L", "K", K_val, value_max
        )

    def update_kl_plane_plot(self):
        if not self.calculator.is_initialized:
            return
        H_val = int(self.kl_ctrl_group["spin"].value())
        hkl_list = self._generate_plane_points("K", "L", "H", H_val)
        results = self.calculator.calculate_structure_factors(hkl_list)
        # compute |F(H_val,0,0)| as color scale max reference
        ref = self.calculator.calculate_structure_factors([[0, 0, 0]])
        value_max = float(np.abs(ref[0])) if len(ref) > 0 else None
        arr = np.array(hkl_list)
        self.kl_plane_visualizer.visualize_plane(
            arr[:, 1], arr[:, 2], np.abs(results), "K", "L", "H", H_val, value_max
        )

    @pyqtSlot()
    def initialize_calculator_inuse(self):
        try:
            params = self.main_window.get_parameters() if self.main_window else None
            cif_path = params.get("cif_file") if params else None
            if not cif_path:
                QMessageBox.warning(
                    self,
                    "Missing CIF",
                    "Please load a valid CIF file in the initialization window first.",
                )
                return
            energy = self.energy_input_inuse.value()
            self.calculator.initialize(cif_path, energy)

            self.status_label_inuse.setText(
                "Status: Calculator initialized successfully"
            )
            self.status_label_inuse.setStyleSheet("color: green; font-weight: bold;")

            # Populate 3D using full HKL cube ranging 0..5
            hkl_list = self._generate_hkl_cube(5)
            results = self.calculator.calculate_structure_factors(hkl_list)
            self.inuse_visualizer3d.visualize_structure_factors(
                hkl_list, np.abs(results)
            )

            # Initialize 2D slices
            self.update_hk_plane_plot()
            self.update_hl_plane_plot()
            self.update_kl_plane_plot()

            # Initialize 3D planes based on current control values
            self.inuse_visualizer3d.set_plane_values(
                L=int(self.hk_ctrl_group["spin"].value()),
                K=int(self.hl_ctrl_group["spin"].value()),
                H=int(self.kl_ctrl_group["spin"].value()),
            )
        except Exception as e:
            QMessageBox.critical(
                self, "Error", f"Failed to initialize calculator: {str(e)}"
            )
            self.status_label_inuse.setText(f"Status: Initialization failed - {str(e)}")
            self.status_label_inuse.setStyleSheet("color: red; font-weight: bold;")

    def _create_customized_tab(self):
        """Create the customized plane subtab with vector inputs and 2D/3D plots."""
        customized_tab = QWidget()
        main_layout = QGridLayout(customized_tab)

        # Left panel: configuration
        left_group = QGroupBox("Configuration")
        left_layout = QFormLayout(left_group)

        # Energy input (in eV)
        self.energy_input = QDoubleSpinBox()
        self.energy_input.setRange(1.0, 100000.0)
        self.energy_input.setDecimals(1)
        self.energy_input.setValue(10000.0)
        self.energy_input.setSuffix(" eV")
        left_layout.addRow("X-ray Energy:", self.energy_input)

        # Compact U, V, Center in one row using text inputs like 110, 010, 000
        uvc_row = QWidget()
        uvc_layout = QHBoxLayout(uvc_row)
        uvc_layout.setContentsMargins(0, 0, 0, 0)
        self.u_line = QLineEdit()
        self.u_line.setPlaceholderText("110")
        self.u_line.setText("110")
        self.v_line = QLineEdit()
        self.v_line.setPlaceholderText("001")
        self.v_line.setText("001")
        self.c_line = QLineEdit()
        self.c_line.setPlaceholderText("000")
        self.c_line.setText("000")
        uvc_layout.addWidget(QLabel("U"))
        uvc_layout.addWidget(self.u_line)
        uvc_layout.addWidget(QLabel("V"))
        uvc_layout.addWidget(self.v_line)
        uvc_layout.addWidget(QLabel("Center"))
        uvc_layout.addWidget(self.c_line)
        left_layout.addRow("", uvc_row)

        # u,v range controls on the same row
        ranges_row = QWidget()
        ranges_layout = QHBoxLayout(ranges_row)
        ranges_layout.setContentsMargins(0, 0, 0, 0)
        self.u_range_spin = QSpinBox()
        self.u_range_spin.setRange(0, 5)
        self.u_range_spin.setValue(3)
        self.v_range_spin = QSpinBox()
        self.v_range_spin.setRange(0, 5)
        self.v_range_spin.setValue(3)
        ranges_layout.addWidget(QLabel("u range"))
        ranges_layout.addWidget(self.u_range_spin)
        ranges_layout.addWidget(QLabel("v range"))
        ranges_layout.addWidget(self.v_range_spin)
        left_layout.addRow("", ranges_row)

        # Initialize button and status
        self.init_btn = QPushButton("Initialize Calculator")
        self.init_btn.clicked.connect(self.initialize_calculator)

        self.status_label = QLabel(
            "Status: Provide CIF in initialization window, then initialize"
        )
        self.status_label.setStyleSheet("color: orange; font-weight: bold;")
        left_layout.addRow("", self.status_label)

        # Update button (kept for manual refresh if needed)
        self.update_plane_btn = QPushButton("Update Plane & Plots")
        self.update_plane_btn.clicked.connect(self.update_customized_plots)
        # Put init and update button on the same row
        buttons_row = QWidget()
        buttons_layout = QHBoxLayout(buttons_row)
        buttons_layout.setContentsMargins(0, 0, 0, 0)
        buttons_layout.addWidget(self.init_btn)
        buttons_layout.addWidget(self.update_plane_btn)
        left_layout.addRow("", buttons_row)

        # React to value changes
        # React to value changes for new compact inputs
        for w in [self.u_range_spin, self.v_range_spin]:
            w.valueChanged.connect(self.update_customized_plots)
        for w in [self.u_line, self.v_line, self.c_line]:
            w.textChanged.connect(self.update_customized_plots)

        # Layout: config (top-left), 2D (bottom-left), 3D (right spanning both rows)
        main_layout.addWidget(left_group, 0, 0, 1, 1)

        # 2D plane visualizer below configuration for more space
        self.custom_plane_visualizer2d = StructureFactorVisualizer2D()
        main_layout.addWidget(self.custom_plane_visualizer2d, 1, 0, 1, 1)

        # Right panel: 3D only
        right_group = QGroupBox("3D Visualization")
        right_layout = QGridLayout(right_group)
        right_layout.addWidget(self.visualizer3d, 0, 0)
        main_layout.addWidget(right_group, 0, 1, 2, 1)

        self.tab_widget.addTab(customized_tab, "Customized plane")

    @pyqtSlot()
    def initialize_calculator(self):
        """Initialize the structure factor calculator for customized plane tab."""
        try:
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
            self.calculator.initialize(cif_path, energy)
            self._calculator_initialized = True
            self.status_label.setText("Status: Calculator initialized successfully")
            self.status_label.setStyleSheet("color: green; font-weight: bold;")
            # Populate initial plots
            self.update_customized_plots()
        except Exception as e:
            QMessageBox.critical(
                self, "Error", f"Failed to initialize calculator: {str(e)}"
            )
            self.status_label.setText(f"Status: Initialization failed - {str(e)}")
            self.status_label.setStyleSheet("color: red; font-weight: bold;")

    def _get_custom_vectors(self):
        """Return U, V, and Center vectors parsed from text inputs like '110'."""

        def parse_hkl(text: str, default: tuple) -> tuple:
            try:
                t = text.strip().replace(" ", "")
                if len(t) != 3 or not all(ch in "0123456789-" for ch in t):
                    return default
                # Support simple forms like 110, 0-11, -10-1 etc.
                # Split by position: allow leading '-' for each digit
                vals = []
                i = 0
                while i < len(t):
                    sign = 1
                    if t[i] == "-":
                        sign = -1
                        i += 1
                    if i >= len(t) or not t[i].isdigit():
                        return default
                    vals.append(sign * int(t[i]))
                    i += 1
                if len(vals) != 3:
                    return default
                return tuple(vals)
            except Exception:
                return default

        U = parse_hkl(self.u_line.text(), (1, 1, 0))
        V = parse_hkl(self.v_line.text(), (0, 0, 1))
        C = parse_hkl(self.c_line.text(), (0, 0, 0))
        return U, V, C

    def update_customized_plots(self):
        """Update 3D scatter (all HKL) with a custom plane overlay and 2D uv plot."""
        try:
            # Always update the plane overlay for immediate feedback
            U, V, C = self._get_custom_vectors()
            u_max = int(self.u_range_spin.value())
            v_max = int(self.v_range_spin.value())
            # symmetric parameter ranges around 0; apply center offset in HKL
            u_min_param = -(u_max // 2)
            u_max_param = u_max - (u_max // 2)
            v_min_param = -(v_max // 2)
            v_max_param = v_max - (v_max // 2)
            try:
                self.visualizer3d.set_custom_plane(
                    U,
                    V,
                    u_min_param,
                    u_max_param,
                    v_min_param,
                    v_max_param,
                    steps=2,
                    center=C,
                )
            except Exception:
                pass

            if not self.calculator.is_initialized:
                return

            # 3D: plot all HKL points 0..5
            hkl_list = self._generate_hkl_cube(5)
            sf_values = self.calculator.calculate_structure_factors(hkl_list)
            self.visualizer3d.visualize_structure_factors(hkl_list, np.abs(sf_values))
            # Re-apply plane overlay after replot
            self.visualizer3d.set_custom_plane(
                U,
                V,
                u_min_param,
                u_max_param,
                v_min_param,
                v_max_param,
                steps=2,
                center=C,
            )

            # 2D: points on the plane using integer combinations of U and V in ranges
            uv_points = []
            hkl_points = []
            # symmetric parameter ranges around 0 with given max steps, shifted by center
            for u in range(u_min_param, u_max_param + 1):
                for v in range(v_min_param, v_max_param + 1):
                    H = C[0] + u * U[0] + v * V[0]
                    K = C[1] + u * U[1] + v * V[1]
                    L = C[2] + u * U[2] + v * V[2]
                    uv_points.append({"u": u, "v": v, "H": H, "K": K, "L": L})
                    hkl_points.append([H, K, L])

            if len(hkl_points) > 0:
                sf_plane = self.calculator.calculate_structure_factors(hkl_points)
                # Reference value for sizing: use |F(0,0,0)| for consistency
                ref = self.calculator.calculate_structure_factors([[0, 0, 0]])
                value_max = (
                    float(np.abs(ref[0]))
                    if len(ref) > 0
                    else (
                        float(np.max(np.abs(sf_plane))) if len(sf_plane) > 0 else None
                    )
                )
                u_label = f"[{U[0]} {U[1]} {U[2]}]"
                v_label = f"[{V[0]} {V[1]} {V[2]}]"
                self.custom_plane_visualizer2d.visualize_uv_plane_points(
                    uv_points,
                    np.abs(sf_plane),
                    u_label,
                    v_label,
                    vector_center=C,
                    value_max=value_max,
                )
        except Exception as e:
            # Keep UI responsive
            print(f"Error updating customized plots: {e}")

    # Legacy custom HKL parser kept for reference (unused in new UI)
    def _parse_custom_hkl(self):
        return []

    def update_results_table(self, hkl_list, sf_values):
        return

    def get_module_instance(self):
        """Get the backend module instance."""
        return self.calculator

    def clear(self):
        """Clear all inputs and results."""
        # Clear customized tab-specific widgets
        self.energy_input.setValue(10000.0)
        self._calculator_initialized = False
        self.status_label.setText(
            "Status: Provide CIF in initialization window, then initialize"
        )
        self.status_label.setStyleSheet("color: orange; font-weight: bold;")

        # Clear visualization
        self.visualizer3d.clear_plot()
        self.custom_plane_visualizer2d.clear_plot()

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
