#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from .visualization import BrillouinVisualizer

# This module would normally import functionality for crystallographic
# calculations using libraries like numpy, scipy, etc.


class BrillouinCalculator:
    """Interface for the Brillouin zone calculator.

    This class handles all the calculations required for the Brillouin zone
    calculator tab. It's a pure Python implementation without PyQt dependencies.
    """

    def __init__(self):
        """Initialize the calculator."""
        self._initialized = False

        # Physical constants
        self.hPlanck = 6.62607015e-34  # Planck's constant [J·s]
        self.c_light = 299792458  # Speed of light [m/s]
        self.e = 1.602176634e-19  # Elementary charge [C]

        # Default lattice parameters
        self.a = 5.0
        self.b = 5.0
        self.c = 13
        self.alpha = 90.0
        self.beta = 90.0
        self.gamma = 90.0
        self.lattice_type = "cubic"
        self.a_star = 2 * np.pi / self.a
        self.b_star = 2 * np.pi / self.b
        self.c_star = 2 * np.pi / self.c

        # X-ray energy and derived quantities
        self.energy = 930  # eV
        self.lambda_A = None  # wavelength in Angstroms
        self.k_in = None  # wavevector magnitude

        # Reciprocal lattice vectors (calculated during initialization)
        self.reciprocal_lattice = None
        self.visualizer = None

    def initialize(
        self, a, b, c, alpha, beta, gamma, energy, lattice_type: str = "cubic"
    ):
        """Initialize with lattice parameters.

        Args:
            a, b, c (float): Lattice constants in Angstroms
            alpha, beta, gamma (float): Lattice angles in degrees
            energy (float): X-ray energy in eV
            lattice_type (str): Lattice type, one of "cubic"... (MARKED: TO BE EXTENDED)
        Returns:
            bool: True if initialization was successful
        """
        try:
            # Store parameters
            self.a = a
            self.b = b
            self.c = c
            self.alpha = alpha
            self.beta = beta
            self.gamma = gamma
            self.energy = energy
            self.lattice_type = lattice_type

            # Calculate wavelength and wavevector
            self.lambda_A = (self.hPlanck * self.c_light) / (energy * self.e) * 1e10
            self.k_in = 1.0 / self.lambda_A

            # Calculate reciprocal lattice
            self._calculate_reciprocal_lattice()

            # Initialize visualizer
            self.visualizer = BrillouinVisualizer(self.reciprocal_lattice)

            self._initialized = True
            return True
        except Exception as e:
            print(f"Error initializing calculator: {str(e)}")
            return False

    def initialize_from_cif(self, cif_file_path, energy):
        """Initialize from a CIF file.
        PROBABLY NOT NEEDED, CONSIDER TO REMOVE IT

        Args:
            cif_file_path (str): Path to the CIF file
            energy (float): X-ray energy in eV

        Returns:
            bool: True if initialization was successful
        """
        try:
            # Simulate loading from CIF file
            # In a real implementation, this would parse the CIF file
            # and extract the lattice parameters

            # For this template, just use default values as if loaded from CIF
            self.a = 5.43  # Silicon lattice constant
            self.b = 5.43
            self.c = 5.43
            self.alpha = 90.0
            self.beta = 90.0
            self.gamma = 90.0
            self.energy = energy

            # Calculate reciprocal lattice
            self._calculate_reciprocal_lattice()

            # Initialize visualizer
            self.visualizer = BrillouinVisualizer(self.reciprocal_lattice)

            self._initialized = True
            return True
        except Exception as e:
            print(f"Error initializing from CIF file: {str(e)}")
            return False

    def _calculate_reciprocal_lattice(self):
        """Calculate the reciprocal lattice vectors."""
        lattice_type = self.lattice_type
        if lattice_type == "cubic":
            self.a_star = 2 * np.pi / self.a
            self.b_star = 2 * np.pi / self.b
            self.c_star = 2 * np.pi / self.c

        # In a real implementation, we would store the full reciprocal lattice vectors
        self.reciprocal_lattice = {
            "a_star": self.a_star,
            "b_star": self.b_star,
            "c_star": self.c_star,
        }

    def calculate_hkl(self, tth, theta, phi, chi):
        """Calculate HKL indices from scattering angles.

        Args:
            tth (float): Scattering angle in degrees
            theta (float): Sample theta rotation in degrees
            phi (float): Sample phi rotation in degrees
            chi (float): Sample chi rotation in degrees

        Returns:
            dict: Dictionary containing H, K, L indices and Q magnitude
        """
        if not self.is_initialized():
            raise ValueError("Calculator not initialized")

        # Calculate Q magnitude
        Q = 2.0 * self.k_in * np.sin(np.radians(tth / 2.0))

        # Calculate delta = theta - (tth/2)
        delta_deg = theta - (tth / 2.0)

        # Build rotation matrices for chi and phi
        chi_rad = np.radians(chi)
        phi_rad = np.radians(phi)

        chi_mat = np.array(
            [
                [1, 0, 0],
                [0, np.cos(chi_rad), -np.sin(chi_rad)],
                [0, np.sin(chi_rad), np.cos(chi_rad)],
            ]
        )

        phi_mat = np.array(
            [
                [np.cos(phi_rad), -np.sin(phi_rad), 0],
                [np.sin(phi_rad), np.cos(phi_rad), 0],
                [0, 0, 1],
            ]
        )

        # Q vector in "theta-2theta coordinates"
        sin_delta = np.sin(np.radians(delta_deg))
        cos_delta = np.cos(np.radians(delta_deg))
        Q_th2th = np.array([Q * sin_delta, 0.0, -Q * cos_delta])

        # Transform to sample coordinates
        q_sample = phi_mat @ chi_mat @ Q_th2th

        # Calculate hkl indices
        h = self.a * q_sample[0]
        k = self.b * q_sample[1]
        l = self.c * q_sample[2]

        Q_rlu = np.sqrt(h**2 + k**2 + l**2)
        return {
            "h": h,
            "k": k,
            "l": l,
            "q": Q_rlu,
            "tth": tth,
            "theta": theta,
            "phi": phi,
            "chi": chi,
        }

    def calculate_angles(self, h, k, l, tth_max=155):
        """Calculate scattering angles from HKL indices.

        CURRENTLY THE CHI IS FIXED TO 0, TO BE EXTENDED

        Args:
            h, k, l (float): HKL indices
            tth_max (float): Maximum allowed scattering angle in degrees

        Returns:
            dict: Dictionary containing scattering angles and minimum energy
        """
        if not self.is_initialized():
            raise ValueError("Calculator not initialized")

        # Convert h,k,l to numpy arrays if they aren't already
        k_h = h / self.a  # k_h has a unit of 1/Å
        k_k = k / self.b  # k_k has a unit of 1/Å
        k_l = l / self.c  # k_l has a unit of 1/Å

        Q_magnitude = np.sqrt(k_h**2 + k_k**2 + k_l**2)
        if Q_magnitude > 2 * self.k_in:
            return {
                "success": False,
                "error": "Energy too low for this momentum transfer. Try to tune down the momentum transfer.",
            }

        tth_rad = 2 * np.arcsin(Q_magnitude / (2 * self.k_in))
        tth = np.degrees(tth_rad)
        if tth > tth_max:
            return {
                "success": False,
                "error": f"Scattering angle (tth) is beyond the maximum allowed value {tth_max}°. Try to tune down the momentum transfer.",
            }

        # calculate Q parallel
        Q_parallel = np.sqrt(k_h**2 + k_k**2)
        theta_rad = tth_rad / 2 - np.arcsin(Q_parallel / Q_magnitude)
        theta = np.degrees(theta_rad)

        # calculate phi
        phi_rad = np.arctan2(k_k, k_h)
        phi = np.degrees(phi_rad)

        return {
            "tth": tth,
            "theta": theta,
            "phi": phi,
            "chi": 0,
            "h": h,
            "k": k,
            "l": l,
            "success": True,
            "error": None,
        }

    def is_initialized(self):
        """Check if the calculator is initialized.

        Returns:
            bool: True if the calculator is initialized
        """
        return self._initialized

    def get_lattice_parameters(self):
        """Get the current lattice parameters.

        Returns:
            dict: Dictionary containing lattice parameters
        """
        return {
            "a": self.a,
            "b": self.b,
            "c": self.c,
            "alpha": self.alpha,
            "beta": self.beta,
            "gamma": self.gamma,
        }
