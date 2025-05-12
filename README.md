<strong style="color: red;"> See [This
site](https://hongxunyang2.github.io/reports/rixs_pre_toolbox/rixs_pre_toolbox) for the REPORTS and the TODO of this project.</strong>


# RIXS Preparation Toolbox

A PyQt5-based application for X-ray spectroscopy preparation, featuring various tools for calculating and visualizing Brillouin zones, scattering geometries, and more.

## Features

- **Global Lattice Parameters**: Initialize crystal structure parameters once and use them across all tools
- **CIF File Support**: Import crystal structure from CIF files or enter parameters manually
- **Brillouin Zone Calculator**: Calculate and visualize Brillouin zones
- **Scattering Geometry Tools**: Convert between angles and HKL indices
- **Modern UI**: Clean, intuitive interface with drag-and-drop support
- **Object-Oriented Architecture**: Hierarchical system of Lab, Sample, and Lattice classes for precise coordinate transformations

## Installation

1. Clone the repository
2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

1. Launch the application
2. Initialize the crystal structure:
   - Enter lattice parameters manually (a, b, c, α, β, γ)
   - Or drag and drop a CIF file
   - Set the X-ray energy
   - Optionally configure sample orientation (roll, pitch, yaw)
3. Use the various tools in the tabs:
   - Brillouin Zone Calculator
   - Angles to HKL converter
   - HKL to Angles converter

## Development

The application is built with:
- Python 3.8+
- PyQt5
- Matplotlib
- NumPy

## Sample Coordinate System Architecture

The application uses a three-level coordinate system architecture:

### 1. Lattice Class
- Manages crystal lattice parameters (a, b, c, α, β, γ)
- Handles real and reciprocal space vectors in the crystal's frame
- Provides methods for unit cell calculations

### 2. Sample Class
- Wraps a Lattice object
- Adds orientation parameters (roll, pitch, yaw) relative to the sample surface
- Transforms vectors between lattice frame and sample frame

### 3. Lab Class
- Wraps a Sample object
- Adds goniometer angles (theta, phi, chi) for sample positioning
- Transforms vectors between sample frame and laboratory frame
- Provides methods for scattering geometry calculations

This hierarchical approach ensures proper coordinate transformations between real/reciprocal space and different reference frames, critical for accurate X-ray scattering calculations.

## License

MIT License

## Module-Tab Communication Architecture

This application uses a clean separation between backend calculations and frontend UI components. Here's how the communication works using the Brillouin Calculator as an example:

### 1. Backend Module (`packages/brillouin_calculator/interface.py`)

The backend module contains a pure Python implementation with no PyQt dependencies:

```python
class BrillouinCalculator:
    def __init__(self):
        # Initialize state
        self._initialized = False
        self.lab = Lab()  # Create a Lab object
        
    def initialize(self, params: dict):
        # Extract parameters from dict
        a, b, c = params.get("a"), params.get("b"), params.get("c")
        alpha, beta, gamma = params.get("alpha"), params.get("beta"), params.get("gamma")
        roll, pitch, yaw = params.get("roll"), params.get("pitch"), params.get("yaw")
        theta, phi, chi = 0.0, 0.0, 0.0  # Default goniometer position
        
        # Initialize the Lab object with parameters
        self.lab.initialize(a, b, c, alpha, beta, gamma, roll, pitch, yaw, theta, phi, chi)
        
    def calculate_hkl(self, tth, theta, phi, chi):
        # Use lab object for calculations
        # Returns result dictionary
```

### 2. Tab Implementation (`packages/gui/tabs/brillouincalculatortab.py`)

The tab implementation creates a UI for the backend module:

```python
class BrillouinCalculatorTab(TabInterface):
    def __init__(self):
        # Create backend instance
        self.calculator = BrillouinCalculator()
        
    def set_parameters(self, params: dict):
        """Set parameters from global settings."""
        self.calculator.initialize(params=params)
```

### 3. Tab Loading (`packages/gui/main_window.py`)

The main window loads all tabs from the registry and passes parameters to each tab.

### Key Communication Patterns

1. **Initialization**: Main window passes parameters to all tabs
2. **Coordinate Transformations**: Lab handles transforms between lab, sample, and lattice frames
3. **Results → UI Update**: Backend returns results, tab updates UI components
4. **Visualization**: Backend generates figures, tab displays them in canvas widgets

## Project Structure

- `packages/`: Contains all functional modules and GUI components
  - `classes/`: Core object model
    - `lab.py`: Lab class for laboratory frame calculations
    - `sample.py`: Sample class for sample frame calculations
    - `lattice.py`: Lattice class for crystal lattice calculations
  - `gui/`: All PyQt5-related code
    - `main_window.py`: Main application window
    - `tabs/`: Tab implementations for each module
  - `visualizer/`: Visualization components
    - `scattering_visualizer.py`: Visualizes scattering geometry
    - `coordinate_visualizer.py`: Visualizes coordinate systems
  - `module_name/`: Each module in its own package
    - `interface.py`: Backend implementation of module functionality
- `config/`: Configuration files for the application
  - `app_config.json`: Main application configuration
  - `tab_registry.json`: Registry of all available tabs
- `static/`: CSS/QSS styling files for the GUI
- `figures/`: Generated figures and visualizations
- `data/`: Data files including example crystal structures

## Dependencies

- Python 3.8+
- PyQt5
- Matplotlib
- NumPy
- CifFile (for crystal structure files)
