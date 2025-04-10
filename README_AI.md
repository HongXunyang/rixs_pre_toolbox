# RIXS Preparation Toolbox - AI Documentation

This document serves as a reference for the AI assistant to understand the structure and design of this project.

## Project Architecture

The project follows a modular architecture to facilitate easy extension and maintenance:

### Core Components

1. **Main Application (`main.py`)**: Entry point that initializes the application and main window
2. **Main Window (`packages/gui/main_window.py`)**: Contains the tab-based interface and loads module GUIs
3. **Module System**: Each functional module is isolated in its own package with a well-defined interface

### Module Design Pattern

Each module follows this pattern:
- **Module Backend (`packages/module_name/`)**: Pure Python implementation without PyQt dependencies
- **Module Interface (`packages/module_name/interface.py`)**: Defines APIs for the GUI to call
- **Module Tab (`packages/gui/tabs/module_name_tab.py`)**: PyQt implementation for the module's GUI

### Tab Registration System

The application uses a tab registry system to dynamically load tabs:
1. Tabs are registered in `config/tab_registry.json`
2. The main window loads tabs dynamically at startup
3. New tabs can be added without modifying existing code

## Directory Structure Details

- `packages/`
  - `gui/`: All PyQt-related code lives here
    - `main_window.py`: Main application window
    - `tabs/`: Individual tab implementations
    - `components/`: Reusable GUI components
    - `utils/`: GUI utility functions
  - `brillouin_calculator/`: Brillouin zone calculator module
  - `brillouin_visualizer/`: Brillouin zone visualizer module
  - `trajectory_planner/`: Trajectory planner module
  
- `config/`
  - `app_config.json`: Global application configuration
  - `tab_registry.json`: Registry of available tabs
  - `module_configs/`: Module-specific configuration files
  
- `static/`
  - `styles.qss`: Main stylesheet for the application
  - `icons/`: Icon resources
  
## Module Extension Process

To add a new module:
1. Create a new package in `packages/` with module implementation
2. Create a tab implementation in `packages/gui/tabs/`
3. Register the tab in `config/tab_registry.json`
4. No changes to existing core code required

## Key Interfaces

Here are the key interfaces between components:
1. Module interface classes in each module package define how the GUI interacts with backend code
2. TabInterface abstract class that all tab implementations must conform to
3. Tab registration mechanism that defines how tabs are dynamically loaded
