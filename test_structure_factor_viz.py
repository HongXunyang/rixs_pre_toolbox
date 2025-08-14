#!/usr/bin/env python3
"""Test script for structure factor visualization"""

import sys
import numpy as np
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget
from packages.visualizer.structure_factor_visualizer import StructureFactorVisualizer

# Test data
test_hkl_list = [
    [1, 0, 0],
    [0, 1, 0], 
    [0, 0, 1],
    [1, 1, 0],
    [1, 0, 1],
    [0, 1, 1],
    [1, 1, 1],
    [2, 0, 0],
    [0, 2, 0],
    [0, 0, 2]
]

# Test structure factor values (some with different magnitudes)
test_sf_values = [10.5, 8.2, 12.0, 6.5, 9.1, 7.8, 15.2, 4.3, 3.9, 11.6]

class TestWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Structure Factor Visualizer Test")
        self.setGeometry(100, 100, 900, 700)
        
        # Create central widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout(central_widget)
        
        # Create visualizer
        self.visualizer = StructureFactorVisualizer()
        layout.addWidget(self.visualizer)
        
        # Test the visualization
        self.test_visualization()
    
    def test_visualization(self):
        print("Testing structure factor visualization...")
        print(f"HKL list: {test_hkl_list}")
        print(f"SF values: {test_sf_values}")
        
        success = self.visualizer.visualize_structure_factors(test_hkl_list, test_sf_values)
        print(f"Visualization success: {success}")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = TestWindow()
    window.show()
    sys.exit(app.exec_())
