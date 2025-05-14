#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import numpy as np
from scipy.optimize import fsolve

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from packages.classes import Lab
from packages.utils import angle_to_matrix
from packages.brillouin_calculator.core import (
    _calculate_angles_factory,
    _calculate_angles_tth_fixed,
    _calculate_hkl,
)


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

        # Initialize sample
        self.lab = Lab()
        self.roll = 0.0
        self.pitch = 0.0
        self.yaw = 0.0

        # X-ray energy and derived quantities
        self.energy = 930  # eV
        self.lambda_A = None  # wavelength in Angstroms
        self.k_in = None  # wavevector magnitude

        # Reciprocal lattice vectors (calculated during initialization)
        self.reciprocal_lattice = None
        self.visualizer = None

    def initialize(
        self,
        params: dict,
    ):
        """Initialize with lattice parameters.

        Args:
            a, b, c (float): Lattice constants in Angstroms
            alpha, beta, gamma (float): Lattice angles in degrees
            energy (float): X-ray energy in eV
            yaw, pitch, roll (float): lattice rotation in degrees
            theta, phi, chi (float): sample rotation in degrees

        Returns:
            bool: True if initialization was successful
        """
        roll = params.get("roll", 0.0)
        pitch = params.get("pitch", 0.0)
        yaw = params.get("yaw", 0.0)
        a, b, c = params.get("a", 4), params.get("b", 4), params.get("c", 12)
        alpha, beta, gamma = (
            params.get("alpha", 90.0),
            params.get("beta", 90.0),
            params.get("gamma", 90.0),
        )
        # the default sample rotation position
        theta, phi, chi = 0.0, 0.0, 0.0
        try:
            # Store parameters
            self.energy = params["energy"]

            # Initialize sample
            self.lab.initialize(
                a,
                b,
                c,
                alpha,
                beta,
                gamma,
                roll,
                pitch,
                yaw,
                theta,
                phi,
                chi,
            )
            # Calculate wavelength and wavevector
            self.lambda_A = (
                (self.hPlanck * self.c_light) / (self.energy * self.e) * 1e10
            )
            self.k_in = 2 * np.pi / self.lambda_A

            self._initialized = True
            return True
        except Exception as e:
            print(f"Error initializing calculator: {str(e)}")
            return False

    def _sample_to_lab_conversion(self, a_vec, b_vec, c_vec):
        """Convert vectors from sample coordinate system to lab coordinate system."""
        # For now, just return the same vectors
        # This should be implemented based on the actual coordinate system conversion
        return a_vec, b_vec, c_vec

    def get_k_magnitude(self, tth):
        return 2.0 * self.k_in * np.sin(np.radians(tth / 2.0))

    def calculate_hkl(self, tth, theta, phi, chi):
        """Calculate HKL from scattering angles.

        Args:
            tth (float): Scattering angle in degrees
            theta (float): Sample theta rotation in degrees
            phi (float): Sample phi rotation in degrees
            chi (float): Sample chi rotation in degrees

        """
        if not self.is_initialized():
            raise ValueError("Calculator not initialized")

        a_vec_lab, b_vec_lab, c_vec_lab = self.lab.get_real_space_vectors()
        return _calculate_hkl(
            self.k_in, tth, theta, phi, chi, a_vec_lab, b_vec_lab, c_vec_lab
        )

    def calculate_angles(
        self,
        H,
        K,
        L,
        fixed_angle,
        fixed_angle_name="chi",
    ):
        """Calculate scattering angles from HKL indices.

        CURRENTLY THE CHI IS FIXED TO 0, TO BE EXTENDED

        Args:
            h, k, l (float): HKL indices

        Returns:
            dict: Dictionary containing scattering angles and minimum energy
        """

        if not self.is_initialized():
            raise ValueError("Calculator not initialized")

        calculate_angles = _calculate_angles_factory(fixed_angle_name)
        a, b, c, alpha, beta, gamma = self.lab.get_lattice_parameters()
        roll, pitch, yaw = self.lab.get_lattice_angles()
        tth_result, theta_result, phi_result, chi_result = calculate_angles(
            self.k_in,
            H,
            K,
            L,
            a,
            b,
            c,
            alpha,
            beta,
            gamma,
            roll,
            pitch,
            yaw,
            fixed_angle,
        )
        return {
            "tth": tth_result,
            "theta": theta_result,
            "phi": phi_result,
            "chi": chi_result,
            "H": H,
            "K": K,
            "L": L,
            "success": True,
            "error": None,
        }

    def calculate_angles_tth_fixed(
        self,
        tth,
        H=0.15,
        K=0.1,
        L=None,
        fixed_angle_name="chi",
        fixed_angle=0.0,
    ):
        a, b, c, alpha, beta, gamma = self.lab.get_lattice_parameters()
        roll, pitch, yaw = self.lab.get_lattice_angles()

        tth_result, theta_result, phi_result, chi_result, momentum = (
            _calculate_angles_tth_fixed(
                self.k_in,
                tth,
                a,
                b,
                c,
                alpha,
                beta,
                gamma,
                roll,
                pitch,
                yaw,
                H,
                K,
                L,
                fixed_angle_name,
                fixed_angle,
            )
        )
        H = momentum if H is None else H
        K = momentum if K is None else K
        L = momentum if L is None else L
        result = {
            "tth": tth_result,
            "theta": theta_result,
            "phi": phi_result,
            "chi": chi_result,
            "H": H,
            "K": K,
            "L": L,
            "success": True,
            "error": None,
        }
        return result

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
        a, b, c, alpha, beta, gamma = self.lab.get_lattice_parameters()
        return {
            "a": a,
            "b": b,
            "c": c,
            "alpha": alpha,
            "beta": beta,
            "gamma": gamma,
        }

    def get_real_space_vectors(self):
        """Get the real space vectors.

        Args:
            frame (str): "sample" or "lab"
        """
        return self.lab.get_real_space_vectors()

