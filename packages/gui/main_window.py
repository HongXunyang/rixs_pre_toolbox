#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
)
from PyQt5.QtCore import Qt, QSize
from PyQt5.QtGui import QIcon


class MainWindow(QMainWindow):
    """Main application window with tab-based interface."""

    def __init__(self):
        super().__init__()

        # Load application configuration
        self.config = self._load_config()

        # Setup window properties
        self.setup_window()

        # Create main UI components
        self.setup_ui()

        # Add tabs from registry
        self.load_tabs()

        # Setup status bar
        self.statusBar().showMessage("Ready")

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
        self.layout.setSpacing(10)  # Add some spacing between widgets

        # Create tab widget
        self.tab_widget = QTabWidget()
        self.tab_widget.setTabPosition(QTabWidget.North)
        self.tab_widget.setMovable(True)
        self.tab_widget.setDocumentMode(True)

        # Add tab widget to layout at (0,0)
        self.layout.addWidget(self.tab_widget, 0, 0)

        # Create menu bar
        self.create_menu_bar()

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

    def load_tabs(self):
        """Load tabs from the registry."""
        registry_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 
                                    "config", "tab_registry.json")

        try:
            with open(registry_path, 'r') as f:
                registry = json.load(f)

            # Sort tabs by order
            tabs_info = sorted(registry.get("tabs", []), key=lambda x: x.get("order", 999))

            # Load each enabled tab
            for tab_info in tabs_info:
                if tab_info.get("enabled", True):
                    self.load_tab(tab_info)

            # Set default tab
            default_tab = self.config.get("default_tab")
            for i in range(self.tab_widget.count()):
                if self.tab_widget.tabText(i) == default_tab:
                    self.tab_widget.setCurrentIndex(i)
                    break

        except Exception as e:
            self.statusBar().showMessage(f"Error loading tabs: {str(e)}")
            if self.config.get("debug_mode", False):
                raise e

    def load_tab(self, tab_info):
        """Load a single tab from the registry."""
        try:
            # Get tab info
            module_name = f"packages.{tab_info['module_package']}"
            tab_class_name = tab_info["tab_class"]

            # Try to load the tab module
            try:
                # First try to load from gui.tabs package
                # Convert class name to lowercase and remove 'Tab' suffix for file name
                file_name = tab_class_name.lower().replace("tab", "")
                module = importlib.import_module(f"packages.gui.tabs.{file_name}")
            except ImportError:
                # If not found, try loading directly from the module package
                module = importlib.import_module(f"{module_name}.gui")

            # Get tab class
            tab_class = getattr(module, tab_class_name)

            # Create tab instance
            tab_instance = tab_class()

            # Add tab to widget
            icon_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 
                                    "static", tab_info.get("icon", ""))

            if os.path.exists(icon_path):
                self.tab_widget.addTab(tab_instance, QIcon(icon_path), tab_info["name"])
            else:
                self.tab_widget.addTab(tab_instance, tab_info["name"])

            # Set tooltip
            index = self.tab_widget.count() - 1
            self.tab_widget.setTabToolTip(index, tab_info.get("description", ""))

        except Exception as e:
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

    def closeEvent(self, event):
        """Handle window close event."""
        # Save session if configured
        if self.config.get("save_session_on_exit", True):
            # TODO: Implement session saving
            pass

        event.accept() 
