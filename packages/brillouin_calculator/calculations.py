#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np


class BrillouinCalculator:
    """Class for performing Brillouin zone calculations."""

    def __init__(self, reciprocal_lattice):
        """Initialize the calculator with reciprocal lattice parameters.

        Args:
            reciprocal_lattice (dict): Dictionary containing reciprocal lattice vectors
        """
        self.reciprocal_lattice = reciprocal_lattice
        self.a_star = reciprocal_lattice["a_star"]
        self.b_star = reciprocal_lattice["b_star"]
        self.c_star = reciprocal_lattice["c_star"]

    def calculate_hkl(self, tth, theta, phi, chi, energy):
        """Calculate HKL indices from scattering angles.

        Args:
            tth (float): Scattering angle in degrees
            theta (float): Sample theta rotation in degrees
            phi (float): Sample phi rotation in degrees
            chi (float): Sample chi rotation in degrees
            energy (float): X-ray energy in eV

        Returns:
            dict: Dictionary containing H, K, L indices and Q magnitude
        """
        # Convert angles to radians
        tth_rad = np.radians(tth)
        theta_rad = np.radians(theta)
        phi_rad = np.radians(phi)
        chi_rad = np.radians(chi)

        # Placeholder calculation of scattering vector
        # These are just placeholder formulas for demonstration
        h = np.sin(tth_rad / 2) * np.cos(theta_rad) * np.cos(phi_rad) / self.a_star
        k = np.sin(tth_rad / 2) * np.cos(theta_rad) * np.sin(phi_rad) / self.b_star
        l = np.sin(tth_rad / 2) * np.sin(theta_rad) / self.c_star

        # Calculate Q magnitude (Å⁻¹)
        wavelength = 12398.42 / energy  # eV to Å
        q = 4 * np.pi * np.sin(tth_rad / 2) / wavelength

        return {
            "h": h,
            "k": k,
            "l": l,
            "q": q,
            "tth": tth,
            "theta": theta,
            "phi": phi,
            "chi": chi,
        }

    def calculate_angles(self, h, k, l, tth_max, fixed_angle, fixed_value, energy):
        """Calculate scattering angles from HKL indices.

        Args:
            h, k, l (float): HKL indices
            tth_max (float): Maximum allowed scattering angle in degrees
            fixed_angle (str): Which angle to fix ('phi' or 'chi')
            fixed_value (float): Value of the fixed angle in degrees
            energy (float): X-ray energy in eV

        Returns:
            dict: Dictionary containing scattering angles and minimum energy
        """
        # Calculate Q magnitude (Å⁻¹)
        q_hkl = np.sqrt(
            (h * self.a_star) ** 2 + (k * self.b_star) ** 2 + (l * self.c_star) ** 2
        )

        # Calculate minimum energy needed for this Q (eV)
        wavelength_max = 4 * np.pi * np.sin(np.radians(tth_max / 2)) / q_hkl
        energy_min = 12398.42 / wavelength_max

        # Calculate wavelength at current energy (Å)
        wavelength = 12398.42 / energy

        # Calculate scattering angle (degrees)
        sin_gamma_2 = q_hkl * wavelength / (4 * np.pi)

        # Check if Q is reachable with current energy
        if sin_gamma_2 > 1:
            # Q not reachable with current energy
            tth = tth_max
            theta = 0.0

            if fixed_angle == "phi":
                phi = fixed_value
                chi = 0.0
            else:
                phi = 0.0
                chi = fixed_value
        else:
            # Q is reachable
            tth = 2 * np.degrees(np.arcsin(sin_gamma_2))

            # Placeholder calculation for other angles
            # In a real implementation, this would involve solving
            # a system of equations to find the angles

            if fixed_angle == "phi":
                phi = fixed_value

                # Placeholder calculations
                theta = np.degrees(
                    np.arctan2(
                        l * self.c_star,
                        np.sqrt((h * self.a_star) ** 2 + (k * self.b_star) ** 2),
                    )
                )

                # Placeholder for chi calculation
                if np.abs(h * self.a_star) + np.abs(k * self.b_star) > 1e-6:
                    chi = np.degrees(np.arctan2(k * self.b_star, h * self.a_star)) - phi
                else:
                    chi = 0.0
            else:
                chi = fixed_value

                # Placeholder calculations
                theta = np.degrees(
                    np.arctan2(
                        l * self.c_star,
                        np.sqrt((h * self.a_star) ** 2 + (k * self.b_star) ** 2),
                    )
                )

                # Placeholder for phi calculation
                if np.abs(h * self.a_star) + np.abs(k * self.b_star) > 1e-6:
                    phi = np.degrees(np.arctan2(k * self.b_star, h * self.a_star)) - chi
                else:
                    phi = 0.0

        return {
            "tth": tth,
            "theta": theta,
            "phi": phi,
            "chi": chi,
            "energy_min": energy_min,
            "h": h,
            "k": k,
            "l": l,
        }
