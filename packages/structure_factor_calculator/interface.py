"""This is a module responsible for calculating the structure factor of a given crystal
structure."""


class StructureFactorCalculator:
    """This is a class responsible for calculating the structure factor of a given crystal
    structure."""

    def __init__(self):
        """Initialize the structure factor calculator."""
        self._initialized = False

    def initialize(self):
        """Initialize the structure factor calculator, the input should be a CIF file"""
        self._initialized = True
        pass

    def is_initialized(self):
        """Check if the structure factor calculator is initialized."""
        return self._initialized
