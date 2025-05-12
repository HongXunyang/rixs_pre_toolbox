"""This is a module responsible for calculating the structure factor of a given crystal
structure."""

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from packages.classes import Lab


class StructureFactorCalculator:
    """This is a class responsible for calculating the structure factor of a given crystal
    structure."""

    def __init__(self):
        """Initialize the structure factor calculator."""
        self._initialized = False

        self.lab = Lab()

        self.energy = 930  # eV

    def initialize(self, params: dict):
        """Initialize the structure factor calculator."""
        self.energy = params.get("energy", 930)  # eV
        a, b, c = params["a"], params["b"], params["c"]
        alpha, beta, gamma = params["alpha"], params["beta"], params["gamma"]
        roll, pitch, yaw = params["roll"], params["pitch"], params["yaw"]
        theta, phi, chi = params["theta"], params["phi"], params["chi"]

        self.lab.initialize(
            a, b, c, alpha, beta, gamma, roll, pitch, yaw, theta, phi, chi
        )

        self._initialized = True

    def is_initialized(self):
        """Check if the structure factor calculator is initialized."""
        return self._initialized

    # def get_lattice_parameters(self) -> dict:
    #     """Get the lattice parameters."""
    #     return self.lab.get_lattice_parameters()

    # def get_real_space_vectors(self) -> dict:
    #     """Get the real space vectors."""
    #     return self.lab.get_real_space_vectors()

    # def get_reciprocal_space_vectors(self) -> dict:
    #     """Get the reciprocal space vectors."""
    #     return self.lab.get_reciprocal_space_vectors()

    def calculate_structure_factor(self, hkl: tuple) -> float:
        """Calculate the structure factor for a given set of hkl indices.

        Args:
            hkl: tuple of integers representing the hkl indices. example: (1, 1, 1)

        Returns:
            float: the structure factor.

        Example:
            >>> calculator = StructureFactorCalculator()
            >>> calculator.initialize(params)
            >>> calculator.calculate_structure_factor((1, 1, 1))
            >>> 1.0
        """
        if not self.is_initialized():
            raise ValueError("Structure factor calculator not initialized")

        # Calculate the structure factor
        structure_factor = None
        return structure_factor

    def calculate_structure_factors(self, hkl_list: list[tuple]) -> list[float]:
        """Calculate the structure factors for a list of hkl indices.

        Args:
            hkl_list: list of tuples of integers representing the hkl indices. example: [(1, 1, 1), (2, 2, 2)]

        Returns:
            list[float]: the structure factors.

        Example:
            >>> calculator = StructureFactorCalculator()
            >>> calculator.initialize(params)
            >>> calculator.calculate_structure_factors([(1, 1, 1), (2, 2, 2)])
            >>> [1.0, 2.0]
        """
        if not self.is_initialized():
            raise ValueError("Structure factor calculator not initialized")

        structure_factors = [None, None]
        return structure_factors
