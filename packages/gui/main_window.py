#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=no-name-in-module, import-error
import os
import json
import importlib
import sys
from PyQt5.QtWidgets import (
    QMainWindow,
    QTabWidget,
    QAction,
    QMenuBar,
    QStatusBar,
    QMessageBox,
    QFileDialog,
    QSplitter,
    QGridLayout,
    QWidget,
    QStackedWidget,
)
from PyQt5.QtCore import Qt, QSize
from PyQt5.QtGui import QIcon

from packages.gui.init_window import InitWindow

class MainWindow(QMainWindow):
    """Main application window with tab-based interface."""

    def __init__(self):
        super().__init__()

        # Load application configuration
        self.config = self._load_config()

        # Initialize lattice parameters
        self.lattice_parameters = None

        # Flag to track if tabs are loaded
        self.tabs_loaded = False

        # Setup window properties
        self.setup_window()

        # Create main UI components
        self.setup_ui()

        # Show initialization window first
        self.show_init_window()

    def setup_window(self):
        """Setup window properties."""
        self.setWindowTitle(self.config.get("app_name", "RIXS Preparation Toolbox"))

        # Set window size
        window_size = self.config.get("window_size", {"width": 1200, "height": 800})
        self.resize(window_size.get("width", 1200), window_size.get("height", 800))

        # Set window icon if available
        icon_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 
                                "static", "icons", "app_icon.png")
        if os.path.exists(icon_path):
            self.setWindowIcon(QIcon(icon_path))

    def setup_ui(self):
        """Setup UI components."""
        # Create central widget and layout
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        self.layout = QGridLayout(self.central_widget)
        self.layout.setContentsMargins(5, 5, 5, 5)
        self.layout.setSpacing(10)

        # Create stacked widget for switching between init and main views
        self.stacked_widget = QStackedWidget()
        self.layout.addWidget(self.stacked_widget, 0, 0)

        # Create initialization window
        self.init_window = InitWindow(self)
        self.stacked_widget.addWidget(self.init_window)

        # Create main tab widget
        self.tab_widget = QTabWidget()
        self.tab_widget.setTabPosition(QTabWidget.West)
        self.tab_widget.setMovable(True)
        self.tab_widget.setDocumentMode(True)
        self.stacked_widget.addWidget(self.tab_widget)

        # Create menu bar
        self.create_menu_bar()

        # Create toolbar
        self.create_toolbar()

        # Create status bar
        self.setStatusBar(QStatusBar())

    def create_menu_bar(self):
        """Create application menu bar."""
        # File menu
        file_menu = self.menuBar().addMenu("&File")

        # File -> Open
        open_action = QAction("&Open", self)
        open_action.setShortcut("Ctrl+O")
        open_action.triggered.connect(self.open_file)
        file_menu.addAction(open_action)

        # File -> Save
        save_action = QAction("&Save", self)
        save_action.setShortcut("Ctrl+S")
        save_action.triggered.connect(self.save_file)
        file_menu.addAction(save_action)

        file_menu.addSeparator()

        # File -> Reset Parameters
        self.reset_action = QAction("&Reset Parameters", self)
        self.reset_action.setShortcut("Ctrl+R")
        self.reset_action.triggered.connect(self.reset_parameters)
        file_menu.addAction(self.reset_action)

        file_menu.addSeparator()

        # File -> Exit
        exit_action = QAction("E&xit", self)
        exit_action.setShortcut("Ctrl+Q")
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)

        # Help menu
        help_menu = self.menuBar().addMenu("&Help")

        # Help -> About
        about_action = QAction("&About", self)
        about_action.triggered.connect(self.show_about)
        help_menu.addAction(about_action)

    def create_toolbar(self):
        """Create application toolbar."""
        # Create toolbar
        self.toolbar = self.addToolBar("Main Toolbar")
        self.toolbar.setMovable(False)
        self.toolbar.setIconSize(QSize(24, 24))

        # Add reset parameters action (reuse the same action from menu)
        reset_icon_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
            "static",
            "icons",
            "reset.png",
        )
        self.reset_action.setIcon(QIcon(reset_icon_path))
        self.reset_action.setToolTip("Reset lattice parameters (Ctrl+R)")
        self.toolbar.addAction(self.reset_action)

        # Add separator
        self.toolbar.addSeparator()

        # Add other toolbar actions here if needed

    def show_init_window(self):
        """Show the initialization window."""
        self.stacked_widget.setCurrentWidget(self.init_window)
        self.statusBar().showMessage("Please initialize lattice parameters")

    def show_main_tabs(self):
        """Show the main tabs after initialization."""
        # Load tabs from registry if not already loaded
        if not self.tabs_loaded:
            self.load_tabs()
            self.tabs_loaded = True

        self.stacked_widget.setCurrentWidget(self.tab_widget)
        self.statusBar().showMessage("Ready")

    def set_lattice_parameters(self, params):
        """Store lattice parameters globally."""
        self.lattice_parameters = params
        # Update all tabs with new parameters
        for i in range(self.tab_widget.count()):
            tab = self.tab_widget.widget(i)
            if hasattr(tab, "set_lattice_parameters"):
                tab.set_lattice_parameters(params)

    def get_lattice_parameters(self):
        """Get the current lattice parameters."""
        return self.lattice_parameters

    def load_tabs(self):
        """Load tabs from the registry."""
        registry_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 
                                    "config", "tab_registry.json")
        print(f"Loading tabs from registry: {registry_path}")

        try:
            with open(registry_path, 'r') as f:
                registry = json.load(f)

            # Sort tabs by order
            tabs_info = sorted(registry.get("tabs", []), key=lambda x: x.get("order", 999))
            print(f"Found {len(tabs_info)} tabs in registry")

            # Load each enabled tab
            for tab_info in tabs_info:
                if tab_info.get("enabled", True):
                    print(f"Loading tab: {tab_info['name']}")
                    self.load_tab(tab_info)

        except Exception as e:
            print(f"Error loading tabs: {str(e)}")
            self.statusBar().showMessage(f"Error loading tabs: {str(e)}")
            if self.config.get("debug_mode", False):
                raise e

    def load_tab(self, tab_info):
        """Load a single tab from the registry."""
        try:
            # Get tab info
            module_name = f"packages.{tab_info['module_package']}"
            tab_class_name = tab_info["tab_class"]
            print(f"Attempting to load tab: {module_name}.{tab_class_name}")

            # Try to load the tab module
            try:
                # First try to load from gui.tabs package
                # Convert class name to lowercase and remove 'Tab' suffix for file name
                file_name = tab_class_name.lower()
                print(f"Trying to load from gui.tabs: {file_name}")
                module = importlib.import_module(f"packages.gui.tabs.{file_name}")
            except ImportError as e:
                print(f"Failed to load from gui.tabs: {str(e)}")
                # If not found, try loading directly from the module package
                print(f"Trying to load from module package: {module_name}.gui")
                module = importlib.import_module(f"{module_name}.gui")

            # Get tab class
            tab_class = getattr(module, tab_class_name)
            print(f"Successfully loaded tab class: {tab_class_name}")

            # Create tab instance
            tab_instance = tab_class(main_window=self)
            print(f"Created tab instance: {tab_instance}")

            # Add tab to widget
            icon_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 
                                    "static", tab_info.get("icon", ""))
            print(f"Loading icon from: {icon_path}")

            if os.path.exists(icon_path):
                self.tab_widget.addTab(tab_instance, QIcon(icon_path), tab_info["name"])
            else:
                self.tab_widget.addTab(tab_instance, tab_info["name"])
            print(f"Added tab to widget: {tab_info['name']}")

            # Set tooltip
            index = self.tab_widget.count() - 1
            self.tab_widget.setTabToolTip(index, tab_info.get("description", ""))

        except Exception as e:
            print(f"Error loading tab {tab_info.get('name')}: {str(e)}")
            self.statusBar().showMessage(f"Error loading tab {tab_info.get('name')}: {str(e)}")
            if self.config.get("debug_mode", False):
                raise e

    def _load_config(self):
        """Load application configuration."""
        config_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 
                                  "config", "app_config.json")
        try:
            with open(config_path, 'r') as f:
                return json.load(f)
        except Exception as e:
            print(f"Error loading config: {str(e)}")
            return {
                "app_name": "RIXS Preparation Toolbox",
                "window_size": {"width": 1200, "height": 800}
            }

    def open_file(self):
        """Open a file."""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Open File", "", "All Files (*);;JSON Files (*.json);;CIF Files (*.cif)"
        )

        if file_path:
            # Pass file to current tab if it has an open_file method
            current_widget = self.tab_widget.currentWidget()
            if hasattr(current_widget, 'open_file'):
                current_widget.open_file(file_path)
            else:
                self.statusBar().showMessage(f"Current tab does not support opening files")

    def save_file(self):
        """Save a file."""
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save File", "", "All Files (*);;JSON Files (*.json)"
        )

        if file_path:
            # Pass file to current tab if it has a save_file method
            current_widget = self.tab_widget.currentWidget()
            if hasattr(current_widget, 'save_file'):
                current_widget.save_file(file_path)
            else:
                self.statusBar().showMessage(f"Current tab does not support saving files")

    def show_about(self):
        """Show about dialog."""
        QMessageBox.about(
            self,
            "About RIXS Preparation Toolbox",
            f"""<b>RIXS Preparation Toolbox</b> v{self.config.get('app_version', '0.1.0')}
            <p>A PyQt5-based application for X-ray spectroscopy preparation.</p>
            <p>Developed as a collaboration project.</p>"""
        )

    def reset_parameters(self):
        """Reset lattice parameters and return to initialization window."""
        # Clear current parameters
        self.lattice_parameters = None

        # Show initialization window
        self.show_init_window()

        # Update status bar
        self.statusBar().showMessage("Please initialize lattice parameters")

    def closeEvent(self, event):
        """Handle window close event."""
        # Save session if configured
        if self.config.get("save_session_on_exit", True):
            # TODO: Implement session saving
            pass

        event.accept() 
