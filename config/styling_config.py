"""Styling configuration for the application."""

import json
import os
from PyQt5.QtGui import QColor, QBrush
class StylingConfig:
    """Styling configuration for the application."""

    def __init__(self):
        """Initialize the styling configuration."""
        self.config = self._load_config()
        
    def _load_config(self):
        """Load the styling configuration."""
        with open(os.path.join(os.path.dirname(__file__), "styling_config.json"), "r") as f:
            return json.load(f)
    
    def get_background_color(self, key):
        """Get the background color for a given key."""
        return self.config["background_colors"][key]
    
    def get_background_qcolor(self, key):
        """Get the background color for a given key."""
        return QColor(*self.get_background_color(key))
    
    def get_background_qbrush(self, key):
        """Get the background brush for a given key."""
        return QBrush(self.get_background_qcolor(key))