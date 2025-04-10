# RIXS Preparation Toolbox

A PyQt5-based application for X-ray spectroscopy preparation with a modular, extensible interface.

## Features

- **Brillouin Zone Calculator**: Calculate and visualize H, K, L indices from angles and vice versa
- **Brillouin Zone Visualizer**: Visualize the Brillouin zone for different crystal structures
- **Trajectory Planner**: Plan measurement trajectories for RIXS experiments
- More features can be added as modular tabs

## Setup

1. Clone this repository
2. Install dependencies:
   ```
   pip install -r requirements.txt
   ```
3. Run the application:
   ```
   python main.py
   ```

## Module-Tab Communication Architecture

This application uses a clean separation between backend calculations and frontend UI components. Here's how the communication works using the Brillouin Calculator as an example:

### 1. Backend Module (`packages/brillouin_calculator/interface.py`)

The backend module contains a pure Python implementation with no PyQt dependencies:

```python
class BrillouinCalculator:
    def __init__(self):
        # Initialize state
        self._initialized = False
        
    def initialize(self, a, b, c, alpha, beta, gamma, energy):
        # Set parameters and calculate 
        # Returns True/False for success/failure
        
    def calculate_hkl(self, gamma, theta, phi, chi):
        # Calculate HKL from angles
        # Returns result dictionary
        
    def visualize_brillouin_zone(self):
        # Create and return a matplotlib Figure
```

### 2. Tab Implementation (`packages/gui/tabs/brillouincalculatortab.py`)

The tab implementation creates a UI for the backend module:

```python
class BrillouinCalculatorTab(TabInterface):
    def __init__(self):
        # Create backend instance
        self.calculator = BrillouinCalculator()
        
    def init_ui(self):
        # Create UI components
        
    def initialize_calculator(self):
        # Collect values from UI
        params = {
            'a': self.a_input.value(),
            # ...other parameters
        }
        # Call backend
        success = self.calculator.initialize(**params)
        # Update UI based on result
        
    def calculate_hkl(self):
        # Collect values from UI
        angles = {
            'gamma': self.gamma_input.value(),
            # ...other angles
        }
        # Call backend
        result = self.calculator.calculate_hkl(**angles)
        # Display results in UI
        
    def update_visualization(self):
        # Get figure from backend
        fig = self.calculator.visualize_brillouin_zone()
        # Display in UI canvas
```

### 3. Tab Loading (`packages/gui/main_window.py`)

The main window loads all tabs from the registry:

```python
def load_tab(self, tab_info):
    # Get tab info from registry
    module_name = f"packages.{tab_info['module_package']}"
    tab_class_name = tab_info["tab_class"]
    
    # Import module and instantiate tab
    module = importlib.import_module(f"packages.gui.tabs.{tab_class_name.lower()}")
    tab_class = getattr(module, tab_class_name)
    tab_instance = tab_class()
    
    # Add tab to the UI
    self.tab_widget.addTab(tab_instance, tab_info["name"])
```

### Key Communication Patterns

1. **Initialization**: Tab creates an instance of the backend module
2. **User Input → Calculation**: UI collects input values and calls backend methods
3. **Results → UI Update**: Backend returns results, tab updates UI components
4. **Visualization**: Backend generates figures, tab displays them in canvas widgets
5. **Error Handling**: Backend returns success/failure, tab shows appropriate messages

This pattern maintains a clean separation of concerns:
- Backend modules focus solely on calculations and logic
- Tab implementations focus on user interface and interaction
- The main window handles tab management and application framework

## Adding New Modules

The application is designed to be easily extendable with new functionality. To add a new module:

1. **Create a new package directory** in the `packages/` folder:
   ```
   packages/your_module_name/
   ```

2. **Implement the backend module** in the package:
   - Create an `__init__.py` file
   - Create an `interface.py` file with a main class that implements your module's functionality
   - Keep this code free of PyQt dependencies

3. **Create a tab implementation** in the GUI package:
   - Create a file in `packages/gui/tabs/yourmoduletab.py`
   - Define a class that inherits from `TabInterface`
   - Implement the UI and connect it to your backend module

4. **Register your tab** in the tab registry:
   - Add an entry to `config/tab_registry.json` with your module information

### Example of a Tab Registry Entry

```json
{
    "id": "your_module",
    "name": "Your Module Name",
    "description": "Description of what your module does",
    "module_package": "your_module_name",
    "tab_class": "YourModuleTab",
    "icon": "icons/your_icon.png",
    "enabled": true,
    "order": 5
}
```

## Project Structure

- `packages/`: Contains all functional modules and GUI components
  - `gui/`: All PyQt5-related code
    - `main_window.py`: Main application window
    - `tabs/`: Tab implementations for each module
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
