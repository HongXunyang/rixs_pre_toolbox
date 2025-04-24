import json


# tips class
class Tips:
    def __init__(self):
        self.tips = {}
        with open("config/tips.json", "r") as f:
            self.tips = json.load(f)

    def tip(self, key):
        return self.tips[key]


def set_tip(widget, tip):
    """Set the tooltip and status tip for a widget."""
    widget.setToolTip(tip)
    widget.setStatusTip(tip)
