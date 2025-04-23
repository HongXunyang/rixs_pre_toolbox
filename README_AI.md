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

## Data Flow

1. **Initialization**
   ```
   User Input → InitWindow → MainWindow → All Tabs
   ```

2. **Parameter Updates**
   ```
   MainWindow (lattice_parameters) → Tab.set_lattice_parameters()
   ```

## Key Features

1. **Global State Management**
   - Lattice parameters stored in `MainWindow`
   - Accessible to all tabs via `get_lattice_parameters()`
   - Updates propagated via `set_lattice_parameters()`

2. **Initialization Process**
   - User enters parameters or drops CIF file
   - Parameters validated and stored
   - Main tabs loaded and initialized
   - Visualization updated

3. **Tab Communication**
   - Tabs implement `set_lattice_parameters()`
   - Parameters updated automatically
   - State preserved across tab switches

## Development Guidelines

1. **Adding New Tabs**
   - Inherit from `TabInterface`
   - Implement `set_lattice_parameters()`
   - Register in tab registry

2. **State Management**
   - Use global parameters from main window
   - Implement state saving/loading
   - Handle parameter updates

3. **UI Design**
   - Follow existing layout patterns
   - Use consistent spacing and margins
   - Support drag-and-drop where appropriate

## Error Handling

1. **Initialization Errors**
   - Validate parameters before storing
   - Show user-friendly error messages
   - Prevent invalid state transitions

2. **Tab Errors**
   - Handle missing parameters gracefully
   - Show appropriate error messages
   - Maintain application stability

## Future Improvements

1. **Parameter Validation**
   - Add more robust validation
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
