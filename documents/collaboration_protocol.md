# RIXS Preparation Toolbox - Collaboration Protocol

This document outlines the protocol for collaboration between UI (frontend) and calculation (backend) developers working on the RIXS Preparation Toolbox.

## Project Structure Overview

```
project/
├── packages/
│   ├── gui/               # Frontend (UI) components
│   │   ├── main_window.py
│   │   ├── tabs/          # Tab implementations
│   │   │   └── *.py
│   ├── module_name/       # Backend calculation modules
│   │   ├── __init__.py
│   │   ├── interface.py   # Primary interface class for backend
│   │   ├── visualization.py  # Visualization support (if needed)
│   │   └── core/          # Core calculations
├── documents/             # Project documentation
```

## 1. Core Principles

1. **Separation of Concerns**
   - Backend: Pure Python calculation logic with no UI dependencies
   - Frontend: User interface components that call backend methods
   
2. **Consistent Interfaces**
   - Backend modules expose a clear, well-documented public API
   - Frontend components handle all user interaction and display logic

3. **Clear Ownership**
   - Backend developer is responsible for calculation accuracy and performance
   - Frontend developer is responsible for user experience and UI components

## 2. Naming Conventions

### 2.1 Module Names
- Backend modules: snake_case, descriptive of functionality (e.g., `brillouin_calculator`)
- Frontend tabs: PascalCase with "Tab" suffix (e.g., `BrillouinCalculatorTab`)

### 2.2 Class Names
- Backend: PascalCase, descriptive of functionality (e.g., `BrillouinCalculator`)
- Frontend: PascalCase, matches tab name (e.g., `BrillouinCalculatorTab`)

### 2.3 Method Names
- Both: snake_case, verb-first to indicate action (e.g., `calculate_hkl`, `update_visualization`)
- Frontend event handlers: include event type (e.g., `on_button_click`, `on_value_changed`)

### 2.4 Variable Names
- Both: snake_case, descriptive without abbreviations unless common
- Parameters: match between frontend/backend for consistency
- Class member variables: prefix with `self.` in methods

### 2.5 Constants
- Both: UPPER_SNAKE_CASE (e.g., `MAX_ANGLE`, `DEFAULT_ENERGY`)

## 3. Data Communication Protocol

### 3.1 Parameter Passing

Backend modules should accept parameters as:

1. **Direct method arguments** for primary calculations
2. **Dictionary with named parameters** for complex operations
3. **Class initialization** for setup parameters

Example:
```python
# Backend
def initialize(self, params: dict):
    self.a = params["a"]
    self.b = params["b"]
    # ...

def calculate_hkl(self, tth, theta, phi, chi):
    # Direct parameters for calculation
    # ...
```

### 3.2 Result Returns

Backend methods should return:

1. **Result dictionary** with named values for multiple returns
2. **Success/error indicators** in the returned dictionary
3. **Specific errors** when calculations fail

Example:
```python
def calculate_hkl(self, tth, theta, phi, chi):
    try:
        # Calculations...
        return {
            "H": h_value,
            "K": k_value,
            "L": l_value,
            "success": True,
            "error": None
        }
    except Exception as e:
        return {
            "success": False,
            "error": str(e)
        }
```

### 3.3 Error Handling

1. Backend should:
   - Validate input parameters
   - Catch and handle calculation errors
   - Return clear error messages

2. Frontend should:
   - Validate user input before calling backend
   - Check success status of backend results
   - Display appropriate error messages to users

## 4. Required Functions

### 4.1 Backend Module Requirements

Each backend module should implement:

1. **Initialization**
   ```python
   def __init__(self):
       # Initialize state variables
       
   def initialize(self, params: dict) -> bool:
       # Set module parameters
       # Return success/failure
   ```

2. **Core Calculation Methods**
   ```python
   def calculate_X(self, param1, param2, ...) -> dict:
       # Perform calculation X
       # Return results dictionary with success flag
   ```

3. **Utility Methods**
   ```python
   def is_initialized(self) -> bool:
       # Check if module is properly initialized
       
   def get_parameters(self) -> dict:
       # Return current parameters
   ```
### 4.2 Frontend Tab Requirements

Each frontend tab should implement:

1. **Initialization**
   ```python
   def __init__(self, main_window=None):
       # Create backend instance
       self.module = ModuleClass()
       # Initialize UI
       
   def init_ui(self):
       # Setup UI components
   ```

2. **Parameter Handling**
   ```python
   def set_parameters(self, params: dict):
       # Update parameters from main window
       self.module.initialize(params)
       # Update UI to reflect new parameters
   ```

3. **Event Handlers**
   ```python
   def on_calculate_button_clicked(self):
       # Get values from UI
       # Call backend method
       # Update UI with results
   ```

4. **Visualization Updates**
   ```python
   def update_visualization(self, results):
       # Update visualization based on calculation results
   ```

## 5. Examples

### 5.1 Backend Module Example

```python
class ExampleCalculator:
    def __init__(self):
        self._initialized = False
        # Default values
        self.a = 5.0
        self.b = 5.0
        
    def initialize(self, params: dict) -> bool:
        try:
            # Store parameters
            self.a = params["a"]
            self.b = params["b"]
            self.energy = params["energy"]
            
            # Perform initialization calculations
            self._calculate_derived_values()
            
            self._initialized = True
            return True
        except Exception as e:
            print(f"Error initializing calculator: {str(e)}")
            return False
            
    def _calculate_derived_values(self):
        # Calculate internal values based on parameters
        pass
        
    def calculate_result(self, input1, input2) -> dict:
        if not self.is_initialized():
            return {"success": False, "error": "Calculator not initialized"}
            
        try:
            # Perform calculation
            result = input1 * self.a + input2 * self.b
            
            return {
                "result": result,
                "input1": input1,
                "input2": input2,
                "success": True,
                "error": None
            }
        except Exception as e:
            return {"success": False, "error": str(e)}
            
    def is_initialized(self) -> bool:
        return self._initialized
```

### 5.2 Frontend Tab Example

```python
class ExampleCalculatorTab(TabInterface):
    def __init__(self, main_window=None):
        # Create backend instance
        self.calculator = ExampleCalculator()
        
        # Initialize UI
        super().__init__(main_window)
        self.init_ui()
        
        # Set parameters if available
        params = self.main_window.get_parameters()
        if params:
            self.set_parameters(params)
        
    def init_ui(self):
        # Create layout
        self.layout = QGridLayout(self)
        
        # Create input fields
        self.input1 = QDoubleSpinBox()
        self.input1.setRange(0, 100)
        self.input1.setValue(1.0)
        
        self.input2 = QDoubleSpinBox()
        self.input2.setRange(0, 100)
        self.input2.setValue(2.0)
        
        # Create calculate button
        self.calculate_button = QPushButton("Calculate")
        self.calculate_button.clicked.connect(self.on_calculate_clicked)
        
        # Create result field
        self.result_field = QLineEdit()
        self.result_field.setReadOnly(True)
        
        # Add widgets to layout
        # ...
        
    def set_parameters(self, params: dict):
        # Initialize calculator with parameters
        success = self.calculator.initialize(params)
        if not success:
            QMessageBox.warning(self, "Warning", "Failed to initialize calculator")
            
    def on_calculate_clicked(self):
        # Get input values
        input1 = self.input1.value()
        input2 = self.input2.value()
        
        # Call backend
        result = self.calculator.calculate_result(input1, input2)
        
        # Check for success
        if result["success"]:
            # Update UI
            self.result_field.setText(str(result["result"]))
        else:
            # Show error
            QMessageBox.warning(self, "Error", result["error"])
```

## 6. Visualization Protocol

For modules requiring visualization:


```

### Frontend Visualization Methods

```python
def update_visualization(self, data):
    # Update visualization canvas with new data
    
def setup_visualization_canvas(self):
    # Setup the visualization canvas in the UI
```

## 7. Git Workflow

1. **Branch Naming**
   - Frontend: `frontend/feature-name`
   - Backend: `backend/feature-name`
   
2. **Commit Messages**
   - Start with area: `[Frontend] Add new input form`
   - Start with area: `[Backend] Implement calculation algorithm`

3. **Pull Requests**
   - Detail interfaces affected
   - Note any API changes
   - Tag collaborator for review

## 8. Documentation Requirements

1. **Backend**
   - Document function parameters and return values
   - Explain calculation algorithms
   - Provide examples
   
2. **Frontend**
   - Document event handlers
   - Explain UI flow
   - Note interactions with backend



By following this protocol, frontend and backend developers can work efficiently together with clear interfaces and expectations. 