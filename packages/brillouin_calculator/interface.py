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

    def calculate_angles(self, h, k, l, tth_max, fixed_angle, fixed_value):
        """Calculate scattering angles from HKL indices.

        Args:
            h, k, l (float): HKL indices
            tth_max (float): Maximum allowed scattering angle in degrees
            fixed_angle (str): Which angle to fix ('phi' or 'chi')
            fixed_value (float): Value of the fixed angle in degrees

        Returns:
            dict: Dictionary containing scattering angles and minimum energy
        """
        if not self.is_initialized():
            raise ValueError("Calculator not initialized")

        # Convert h,k,l to numpy arrays if they aren't already
        h = np.atleast_1d(h)
        k = np.atleast_1d(k)
        l = np.atleast_1d(l)

        # Initialize output arrays
        n = len(h)
        tth = np.zeros(n)
        theta = np.zeros(n)
        phi = np.zeros(n)
        chi = np.zeros(n)  # Initialize chi array
        delta = np.zeros(n)
        theta_out = np.zeros(n)
        flag = False  # Flag to indicate if Q is reachable

        # First calculate phi from h and k
        for i in range(n):
            # Calculate phi for each pair of values
            if h[i] != 0:
                phi[i] = np.degrees(np.arctan(np.abs(k[i] / self.b / (h[i] / self.a))))
            else:
                # Handle Gamma point (h=0, k=0)
                if k[i] == 0:
                    # We're at Gamma point, phi is undefined
                    # Use a small value or neighboring point as in MATLAB code
                    if i > 0:
                        haux = h[i - 1] / 1e5
                        kaux = k[i - 1] / 1e5
                    else:
                        haux = h[i + 1] / 1e5 if i < n - 1 else 1e-5
                        kaux = k[i + 1] / 1e5 if i < n - 1 else 1e-5

                    phi[i] = np.degrees(np.arctan(np.abs(kaux / haux)))
                else:
                    # k is non-zero but h is zero
                    phi[i] = 90.0

            # Calculate total momentum transfer
            q_hkl = np.sqrt(
                (h[i] / self.a) ** 2 + (k[i] / self.b) ** 2 + (l / self.c) ** 2
            )
            q_parallel = np.sign(h[i] + 1e-15) * np.sqrt(
                (h[i] / self.a) ** 2 + (k[i] / self.b) ** 2
            )
            q_perp = l / self.c

            # Check if Q is reachable with current energy/max scattering angle
            if q_hkl / (2 * self.k_in * np.sin(np.radians(tth_max / 2))) > 1:
                # Q is not reachable
                flag = True
                continue
            else:
                # Q is reachable, calculate angles
                tth[i] = 2 * np.degrees(np.arcsin(q_hkl / (2 * self.k_in)))
                theta[i] = (
                    np.degrees(np.arctan(q_parallel / np.abs(q_perp))) + tth[i] / 2
                )
                delta[i] = theta[i] - tth[i] / 2
                theta_out[i] = tth[i] - theta[i]

        # Apply fixed angle constraints (phi or chi fixed)
        if fixed_angle == "phi":
            if not flag:
                # Q is reachable, enforce fixed phi value
                phi = np.full_like(phi, fixed_value)
                # Chi needs to be calculated based on fixed phi
                # This is a simplification - in a real implementation,
                # you'd calculate chi based on the fixed phi value
        else:  # fixed_angle == "chi"
            if not flag:
                # Q is reachable, enforce fixed chi value
                chi = np.full_like(phi, fixed_value)

        # Calculate minimum energy needed for this Q
        q_hkl_max = np.max(
            np.sqrt((h / self.a) ** 2 + (k / self.b) ** 2 + (l / self.c) ** 2)
        )
        wavelength_max = 4 * np.pi * np.sin(np.radians(tth_max / 2)) / q_hkl_max
        energy_min = 12398.42 / wavelength_max  # Convert wavelength to energy in eV

        # Return as single values if input was scalar
        if len(h) == 1:
            return {
                "tth": tth[0],
                "theta": theta[0],
                "phi": phi[0],
                "chi": chi[0],
                "delta": delta[0],
                "theta_out": theta_out[0],
                "energy_min": energy_min,
                "h": h[0],
                "k": k[0],
                "l": l,
                "flag": flag,
            }
        else:
            return {
                "tth": tth,
                "theta": theta,
                "phi": phi,
                "chi": chi,
                "delta": delta,
                "theta_out": theta_out,
                "energy_min": energy_min,
                "h": h,
                "k": k,
                "l": l,
                "flag": flag,
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
