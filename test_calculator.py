#!/usr/bin/env python3
"""Test the structure factor calculator to see what data it returns"""

import numpy as np
from packages.structure_factor_calculator.interface import StructureFactorCalculator, DEFAULT_HKL_LIST

# Test with a CIF file
calculator = StructureFactorCalculator()

# Check if the test CIF file exists
import os
cif_file = "data/nacl.cif"
if os.path.exists(cif_file):
    print(f"Using CIF file: {cif_file}")
    
    # Initialize calculator
    calculator.initialize(cif_file, 10000.0)  # 10 keV
    print("Calculator initialized successfully")
    
    # Calculate structure factors
    results = calculator.calculate_structure_factors()
    
    print(f"Default HKL list: {DEFAULT_HKL_LIST}")
    print(f"Results type: {type(results)}")
    print(f"Results shape: {np.array(results).shape}")
    print(f"Results dtype: {np.array(results).dtype}")
    print(f"First few results: {results[:3]}")
    print(f"Results magnitude: {np.abs(results)[:3]}")
    
    # Check if results are complex
    if hasattr(results[0], 'real'):
        print("Results are complex numbers")
        print(f"Real parts: {[r.real for r in results[:3]]}")
        print(f"Imaginary parts: {[r.imag for r in results[:3]]}")
    else:
        print("Results are real numbers")
        
else:
    print(f"CIF file {cif_file} not found")
