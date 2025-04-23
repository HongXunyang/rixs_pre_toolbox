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


class MatplotlibCanvas(FigureCanvas):
    """Matplotlib canvas for embedding in Qt applications."""

    def __init__(self, width=6, height=5, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111, projection="3d")
        super(MatplotlibCanvas, self).__init__(self.fig)
        self.fig.tight_layout()


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

    def __init__(self):
        # Create backend instance
        self.calculator = BrillouinCalculator()

        # Initialize UI
        super().__init__()

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

    def set_lattice_parameters(self, params):
        """Set lattice parameters from global settings."""
        try:
            if params.get("cif_file"):
                # Use CIF file
                success = self.calculator.initialize_from_cif(
                    params["cif_file"], params["energy"]
                )
            else:
                # Use manual parameters
                success = self.calculator.initialize(
                    a=params["a"],
                    b=params["b"],
                    c=params["c"],
                    alpha=params["alpha"],
                    beta=params["beta"],
                    gamma=params["gamma"],
                    energy=params["energy"],
                )

            if not success:
                QMessageBox.warning(self, "Error", "Failed to initialize calculator!")

        except Exception as e:
            QMessageBox.critical(
                self, "Error", f"Error initializing calculator: {str(e)}"
            )

    def create_angles_to_hkl_tab(self):
        """Create tab for angles to HKL calculation."""
        angles_tab = QWidget()
        angles_layout = QVBoxLayout(angles_tab)

        # Input form
        form_group = QGroupBox("Scattering Angles")
        form_layout = QFormLayout(form_group)

        self.gamma_angle_input = QDoubleSpinBox()
        self.gamma_angle_input.setRange(0.0, 180.0)
        self.gamma_angle_input.setValue(30.0)
        self.gamma_angle_input.setSuffix(" °")
        form_layout.addRow("γ (scattering):", self.gamma_angle_input)

        self.theta_angle_input = QDoubleSpinBox()
        self.theta_angle_input.setRange(-180.0, 180.0)
        self.theta_angle_input.setValue(0.0)
        self.theta_angle_input.setSuffix(" °")
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

        angles_layout.addWidget(form_group)

        # Calculate button
        calculate_button = QPushButton("Calculate HKL")
        calculate_button.clicked.connect(self.calculate_hkl)
        angles_layout.addWidget(calculate_button)

        # Results group
        results_group = QGroupBox("Results")
        results_layout = QFormLayout(results_group)

        self.h_result = QLineEdit()
        self.h_result.setReadOnly(True)
        results_layout.addRow("H:", self.h_result)

        self.k_result = QLineEdit()
        self.k_result.setReadOnly(True)
        results_layout.addRow("K:", self.k_result)

        self.l_result = QLineEdit()
        self.l_result.setReadOnly(True)
        results_layout.addRow("L:", self.l_result)

        self.q_result = QLineEdit()
        self.q_result.setReadOnly(True)
        results_layout.addRow("Q (Å⁻¹):", self.q_result)

        angles_layout.addWidget(results_group)

        # Add spacer
        angles_layout.addStretch()

        # Add to tab widget
        self.tab_widget.addTab(angles_tab, "Angles → HKL")

    def create_hkl_to_angles_tab(self):
        """Create tab for HKL to angles calculation."""
        hkl_tab = QWidget()
        hkl_layout = QVBoxLayout(hkl_tab)

        # Input form
        form_group = QGroupBox("HKL Indices")
        form_layout = QFormLayout(form_group)

        self.h_input = QDoubleSpinBox()
        self.h_input.setRange(-10.0, 10.0)
        self.h_input.setDecimals(3)
        self.h_input.setValue(1.0)
        form_layout.addRow("H:", self.h_input)

        self.k_input = QDoubleSpinBox()
        self.k_input.setRange(-10.0, 10.0)
        self.k_input.setDecimals(3)
        self.k_input.setValue(0.0)
        form_layout.addRow("K:", self.k_input)

        self.l_input = QDoubleSpinBox()
        self.l_input.setRange(-10.0, 10.0)
        self.l_input.setDecimals(3)
        self.l_input.setValue(0.0)
        form_layout.addRow("L:", self.l_input)

        hkl_layout.addWidget(form_group)

        # Constraints group
        constraints_group = QGroupBox("Constraints")
        constraints_layout = QFormLayout(constraints_group)

        self.gamma_max_input = QDoubleSpinBox()
        self.gamma_max_input.setRange(0.0, 180.0)
        self.gamma_max_input.setValue(120.0)
        self.gamma_max_input.setSuffix(" °")
        constraints_layout.addRow("γ max:", self.gamma_max_input)

        # Fixed angle selection
        fixed_angle_layout = QHBoxLayout()
        self.fixed_phi_radio = QRadioButton("Fix φ")
        self.fixed_chi_radio = QRadioButton("Fix χ")
        self.fixed_phi_radio.setChecked(True)
        fixed_angle_layout.addWidget(self.fixed_phi_radio)
        fixed_angle_layout.addWidget(self.fixed_chi_radio)
        constraints_layout.addRow("Fixed angle:", fixed_angle_layout)

        # Fixed angle value
        self.fixed_angle_value = QDoubleSpinBox()
        self.fixed_angle_value.setRange(-180.0, 180.0)
        self.fixed_angle_value.setValue(0.0)
        self.fixed_angle_value.setSuffix(" °")
        constraints_layout.addRow("Fixed value:", self.fixed_angle_value)

        hkl_layout.addWidget(constraints_group)

        # Calculate button
        calculate_button = QPushButton("Calculate Angles")
        calculate_button.clicked.connect(self.calculate_angles)
        hkl_layout.addWidget(calculate_button)

        # Results group
        results_group = QGroupBox("Results")
        results_layout = QFormLayout(results_group)

        self.gamma_result = QLineEdit()
        self.gamma_result.setReadOnly(True)
        results_layout.addRow("γ (scattering):", self.gamma_result)

        self.theta_result = QLineEdit()
        self.theta_result.setReadOnly(True)
        results_layout.addRow("θ:", self.theta_result)

        self.phi_result = QLineEdit()
        self.phi_result.setReadOnly(True)
        results_layout.addRow("φ:", self.phi_result)

        self.chi_result = QLineEdit()
        self.chi_result.setReadOnly(True)
        results_layout.addRow("χ:", self.chi_result)

        self.energy_min_result = QLineEdit()
        self.energy_min_result.setReadOnly(True)
        results_layout.addRow("Min. Energy:", self.energy_min_result)

        hkl_layout.addWidget(results_group)

        # Add spacer
        hkl_layout.addStretch()

        # Add to tab widget
        self.tab_widget.addTab(hkl_tab, "HKL → Angles")

    def create_visualization_area(self):
        """Create area for visualization of results."""
        vis_group = QGroupBox("Visualization")
        vis_layout = QVBoxLayout(vis_group)

        # Add matplotlib canvas
        self.canvas = MatplotlibCanvas(width=6, height=5)
        vis_layout.addWidget(self.canvas)

        # Add to main layout
        self.layout.addWidget(vis_group)

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
            result = self.calculator.calculate_hkl(
                gamma=self.gamma_angle_input.value(),
                theta=self.theta_angle_input.value(),
                phi=self.phi_angle_input.value(),
                chi=self.chi_angle_input.value(),
            )

            # Update results
            self.h_result.setText(f"{result['h']:.4f}")
            self.k_result.setText(f"{result['k']:.4f}")
            self.l_result.setText(f"{result['l']:.4f}")
            self.q_result.setText(f"{result['q']:.4f}")

            # Update visualization
            self.update_visualization(result)

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error calculating HKL: {str(e)}")

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

            # Determine which angle is fixed
            fixed_angle = "phi" if self.fixed_phi_radio.isChecked() else "chi"
            fixed_value = self.fixed_angle_value.value()

            # Calculate angles
            result = self.calculator.calculate_angles(
                h=self.h_input.value(),
                k=self.k_input.value(),
                l=self.l_input.value(),
                gamma_max=self.gamma_max_input.value(),
                fixed_angle=fixed_angle,
                fixed_value=fixed_value,
            )

            # Update results
            self.gamma_result.setText(f"{result['gamma']:.4f}°")
            self.theta_result.setText(f"{result['theta']:.4f}°")
            self.phi_result.setText(f"{result['phi']:.4f}°")
            self.chi_result.setText(f"{result['chi']:.4f}°")

            if "energy_min" in result:
                self.energy_min_result.setText(f"{result['energy_min']:.1f} eV")
            else:
                self.energy_min_result.setText("N/A")

            # Update visualization
            self.update_visualization(result)

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error calculating angles: {str(e)}")

    def update_visualization(self, data=None):
        """Update the visualization with current data."""
        # Call the calculator's visualization method
        if data is None:
            # Just draw the basic Brillouin zone
            fig = self.calculator.visualize_brillouin_zone()
        else:
            # Draw with scattering geometry
            fig = self.calculator.visualize_scattering_geometry(data)

        if fig:
            # Replace the current figure in the canvas
            self.canvas.fig.clear()
            for ax in fig.get_axes():
                ax.get_figure().delaxes(ax)
                ax.figure = self.canvas.fig
                self.canvas.fig.add_axes(ax)

            self.canvas.draw()

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
