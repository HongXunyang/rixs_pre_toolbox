# Brillouin Calculator Module Documentation

## Overview

The `brillouin_calculator` module provides:
- Conversion between scattering angles and HKL indices
- Reciprocal lattice calculations


## Class: `BrillouinCalculator`

### Core Methods


#### `initialize(self, params: dict) -> bool`
Sets up the calculator with lattice parameters. This method is required for all the backend classes. 

**Parameters**: Dictionary containing `a`, `b`, `c`, `alpha`, `beta`, `gamma`, `energy`, and optionally `e_H`, `e_K`, `e_L`.

**Returns**: Boolean success indicator.

#### `calculate_hkl(self, tth, theta, phi, chi) -> dict`
Converts scattering angles to HKL indices.

**Parameters**: `tth`, `theta`, `phi`, `chi` angles in degrees.

**Returns**: Dictionary with `H`, `K`, `L`, angles, and success status.

#### `calculate_angles(self, H_crystal, K_crystal, L_crystal, fixed_angle, fixed_angle_name="chi") -> dict`
Calculates scattering angles from HKL with one angle fixed.

**Parameters**: 
- `H_crystal`, `K_crystal`, `L_crystal`: HKL indices
- `fixed_angle`: Value of fixed angle
- `fixed_angle_name`: Which angle to fix ("chi" or "phi")

**Returns**: Dictionary with calculated angles and success status.

#### `calculate_angles_tth_fixed(self, tth, H_crystal=0.15, K_crystal=0.1, L_crystal=None, fixed_angle_name="chi", fixed_angle=0.0) -> dict`
Calculates angles from two HKL values with fixed tth.

**Parameters**:
- `tth`: Fixed scattering angle
- `H_crystal`, `K_crystal`, `L_crystal`: HKL indices (one can be None)
- Fixed angle parameters

**Returns**: Dictionary with calculated values.

#### `is_initialized(self) -> bool`
Checks if calculator is ready for use.

#### `get_lattice_parameters(self) -> dict`
Returns current lattice parameters.


## Usage Example

```python
from packages.brillouin_calculator.interface import BrillouinCalculator

# Create and initialize
calculator = BrillouinCalculator()
params = {
    "a": 5.43, "b": 5.43, "c": 5.43,
    "alpha": 90.0, "beta": 90.0, "gamma": 90.0,
    "energy": 930.0
}
calculator.initialize(params)

# Calculate HKL from angles
result = calculator.calculate_hkl(tth=150.0, theta=75.0, phi=0.0, chi=0.0)
if result["success"]:
    print(f"H={result['H']:.4f}, K={result['K']:.4f}, L={result['L']:.4f}")
```

## Error Handling

All calculation methods return a dictionary with `success` and `error` keys. Always check the `success` key before using results. 