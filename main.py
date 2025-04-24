#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=no-name-in-module, import-error
import sys
import os
from PyQt5.QtWidgets import QApplication
from PyQt5.QtCore import Qt
from packages.gui.main_window import MainWindow


def load_stylesheet():
    """Load QSS stylesheet for the application."""
    stylesheet_path = os.path.join(os.path.dirname(__file__), 'static', 'styles.qss')
    if os.path.exists(stylesheet_path):
        with open(stylesheet_path, 'r') as f:
            return f.read()
    return ""


def main():
    """Main application entry point."""
    # Enable high DPI scaling
    QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)
    QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True)
    
    # Create application
    app = QApplication(sys.argv)
    app.setApplicationName("RIXS Preparation Toolbox")
    
    # Apply stylesheet
    app.setStyleSheet(load_stylesheet())
    
    # Create and show main window
    main_window = MainWindow()
    main_window.show()
    
    # Start application event loop
    sys.exit(app.exec_())


if __name__ == "__main__":
    main() 
