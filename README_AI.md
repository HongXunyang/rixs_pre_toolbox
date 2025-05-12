# RIXS Preparation Toolbox - AI Documentation

## Architecture Overview

The application follows a modular architecture with the following key components:

1. **Main Window (`MainWindow`)**
   - Central application window
   - Manages global state (lattice parameters)
   - Handles tab switching and initialization
   - Uses `QStackedWidget` to switch between initialization and main views

2. **Initialization Window (`InitWindow`)**
   - First window shown to users
   - Handles crystal structure initialization
   - Supports both manual parameter entry and CIF file import
   - Drag-and-drop support for CIF files

3. **Tab Interface (`TabInterface`)**
   - Base class for all tool tabs
   - Provides common functionality
   - Implements state saving/loading

4. **Tool Tabs**
   - Brillouin Zone Calculator
   - Angles to HKL converter
   - HKL to Angles converter
   - Each tab receives lattice parameters from main window

5. **Core Object Model**
   - **Lattice Class**: Manages crystal lattice parameters and vectors
   - **Sample Class**: Adds sample orientation relative to crystal lattice
   - **Lab Class**: Adds goniometer angles for sample positioning in lab

## Data Flow

1. **Initialization**
   ```
   User Input → InitWindow → MainWindow → All Tabs
   ```

2. **Parameter Updates**
   ```
   MainWindow (parameters) → Tab.set_parameters() → Backend.initialize(params)
   ```

3. **Coordinate Transformations**
   ```
   Lattice (crystal frame) → Sample (sample frame) → Lab (lab frame)
   ```

## Coordinate Systems

1. **Crystal (Lattice) Frame**
   - Origin at lattice origin
   - X-axis along a
   - XY-plane contains b
   - Managed by `Lattice` class

2. **Sample Frame**
   - Origin at sample center
   - Z-axis normal to sample surface
   - XY-plane is sample surface
   - Managed by `Sample` class
   - Transformed from crystal frame by roll, pitch, yaw angles

3. **Laboratory Frame**
   - Origin at beam-sample intersection
   - Z-axis along incident beam direction
   - XY-plane perpendicular to beam
   - Managed by `Lab` class
   - Transformed from sample frame by theta, phi, chi angles

## Key Features

1. **Global State Management**
   - Parameters stored in `MainWindow`
   - Passed to tabs via `set_parameters()`
   - Parameters used to initialize backend calculators

2. **Initialization Process**
   - User enters parameters or drops CIF file
   - Parameters validated and stored in MainWindow
   - Main tabs loaded and receive parameters
   - Each tab initializes its backend with parameters

3. **Coordinate Transformations**
   - Hierarchical transformation chain
   - Real and reciprocal space vector calculations
   - Automatic updates when angles change

4. **Scattering Calculations**
   - Lab class handles momentum transfer calculations
   - Proper coordinate transforms for accurate results
   - Consistent calculations across all tools

## Development Guidelines

1. **Adding New Tabs**
   - Inherit from `TabInterface`
   - Implement `set_parameters(params)`
   - Create backend calculator that uses the Lab class

2. **State Management**
   - Use parameters dict from main window
   - Initialize Lab object with parameters
   - Update Lab object when parameters change

3. **UI Design**
   - Follow existing layout patterns
   - Use consistent spacing and margins
   - Support drag-and-drop where appropriate

4. **Coordinate Handling**
   - Always use the appropriate coordinate system
   - Document which frame vectors are in
   - Use Lab, Sample, Lattice objects for transformations

## Error Handling

1. **Initialization Errors**
   - Validate parameters before creating objects
   - Use try/except in initialization methods
   - Provide clear error messages for invalid parameters

2. **Tab Errors**
   - Check for initialized calculator
   - Handle missing parameters gracefully
   - Show appropriate error messages

## Future Improvements

1. **Parameter Validation**
   - Add more robust validation for lattice parameters
   - Support for different crystal systems
   - Unit conversion tools

2. **State Management**
   - Session saving/loading
   - Parameter presets
   - Undo/redo functionality

3. **UI Enhancements**
   - Dark mode support
   - Customizable layouts
   - More visualization options

4. **Coordinate Systems**
   - Additional reference frames
   - More flexible goniometer configurations
   - Support for custom sample environments
