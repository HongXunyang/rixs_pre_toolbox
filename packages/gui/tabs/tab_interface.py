#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from PyQt5.QtWidgets import QWidget, QVBoxLayout


class TabInterface(QWidget):
    """Base class for all tab implementations.
    
    All tab implementations should inherit from this class and implement
    the required methods to ensure consistent behavior across tabs.
    """
    
    def __init__(self):
        super().__init__()
        
        # Create main layout
        self.layout = QVBoxLayout(self)
        self.layout.setContentsMargins(10, 10, 10, 10)
        
        # Initialize UI components
        self.init_ui()
    
    def init_ui(self):
        """Initialize UI components.
        
        This method should be implemented by subclasses to create all UI components.
        """
        raise NotImplementedError("Subclasses must implement init_ui()")
    
    def open_file(self, file_path):
        """Handle opening a file.
        
        Args:
            file_path (str): Path to the file to open
            
        Returns:
            bool: True if file was opened successfully, False otherwise
        """
        # Default implementation - subclasses should override if needed
        return False
    
    def save_file(self, file_path):
        """Handle saving to a file.
        
        Args:
            file_path (str): Path to save the file to
            
        Returns:
            bool: True if file was saved successfully, False otherwise
        """
        # Default implementation - subclasses should override if needed
        return False
    
    def clear(self):
        """Clear all inputs and results."""
        # Default implementation - subclasses should override if needed
        pass
    
    def get_module_instance(self):
        """Get the backend module instance for this tab.
        
        Returns:
            object: Instance of the backend module
        """
        # Default implementation - subclasses should override if needed
        return None

    def get_state(self):
        """Get the current state of the tab for session saving.
        
        Returns:
            dict: Current state of the tab
        """
        # Default implementation - subclasses should override if needed
        return {}
    
    def set_state(self, state):
        """Restore tab state from saved session.
        
        Args:
            state (dict): State to restore
            
        Returns:
            bool: True if state was restored successfully, False otherwise
        """
        # Default implementation - subclasses should override if needed
        return False 