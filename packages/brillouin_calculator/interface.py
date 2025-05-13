#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# import partial
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from functools import partial
import numpy as np
from packages.classes import Lab
from packages.utils import angle_to_matrix

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

        # Calculate momentum transfer magnitude
        k_magnitude = self.get_k_magnitude(tth)
        # Calculate delta = theta - (tth/2)
        delta = -(tth / 2.0)
        sin_delta = np.sin(np.radians(delta))
        cos_delta = np.cos(np.radians(delta))
        # momentum transfer at theta, phi, chi = 0
        k_vec_initial = np.array(
            [k_magnitude * sin_delta, 0.0, -k_magnitude * cos_delta]
        )
        # rotation of the beam is the reverse rotation of the sample, thus the transpose
        rotation_matrix = angle_to_matrix(theta, phi, chi).T
        # momentum transfer at non-zero theta, phi, chi
        k_vec_lab = rotation_matrix @ k_vec_initial
        a_vec_lab, b_vec_lab, c_vec_lab = self.lab.get_real_space_vectors()

        # calculate HKL
        H = np.dot(k_vec_lab, a_vec_lab) / (2 * np.pi)
        K = np.dot(k_vec_lab, b_vec_lab) / (2 * np.pi)
        L = np.dot(k_vec_lab, c_vec_lab) / (2 * np.pi)
        return {
            "H": H,
            "K": K,
            "L": L,
            "tth": tth,
            "theta": theta,
            "phi": phi,
            "chi": chi,
            "success": True,
            "error": None,
        }

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
        a_star_vec_lab, b_star_vec_lab, c_star_vec_lab = (
            self.lab.get_reciprocal_space_vectors()
        )
        return calculate_angles(
            self.k_in,
            H,
            K,
            L,
            a_star_vec_lab,
            b_star_vec_lab,
            c_star_vec_lab,
            fixed_angle,
        )

    def calculate_angles_tth_fixed(
        self,
        tth,
        H=0.15,
        K=0.1,
        L=None,
        fixed_angle_name="chi",
        fixed_angle=0.0,
    ):
        calculate_angles = _calculate_angles_factory("tth")
        a_vec_lab, b_vec_lab, c_vec_lab = self.lab.get_real_space_vectors()
        a_star_vec_lab, b_star_vec_lab, c_star_vec_lab = (
            self.lab.get_reciprocal_space_vectors()
        )
        results = calculate_angles(
            self.k_in,
            tth,
            H,
            K,
            L,
            a_vec_lab,
            b_vec_lab,
            c_vec_lab,
            a_star_vec_lab,
            b_star_vec_lab,
            c_star_vec_lab,
            fixed_angle_name,
            fixed_angle,
        )
        return results

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


def _calculate_angles_factory(fixed_angle_name):
    if fixed_angle_name == "chi":
        return _calculate_angles_chi_fixed
    elif fixed_angle_name == "phi":
        return _calculate_angles_phi_fixed
    elif fixed_angle_name == "tth":
        return _calculate_angles_tth_fixed


def _calculate_angles_chi_fixed(
    k_in, H, K, L, a_star_vec_lab, b_star_vec_lab, c_star_vec_lab, chi
):
    """Calculate scattering angles from HKL, with chi fixed."""

    # convert crystal coordinates to lab coordinates
    epsilon = 1e-10

    # get momentum transfer vector in lab coordinates
    k_vec_lab = H * a_star_vec_lab + K * b_star_vec_lab + L * c_star_vec_lab
    print(f"k_vec_lab: {k_vec_lab}")
    cos_chi = np.cos(np.radians(chi))
    sin_chi = np.sin(np.radians(chi))
    k_magnitude = np.linalg.norm(k_vec_lab)
    tth_rad = 2 * np.arcsin(k_magnitude / (2 * k_in))
    tth = np.degrees(tth_rad)  # tth is done calculating

    # delta == theta - (tth/2) is the angle between the momentum transfer vector and sample surface
    delta_rad = np.sign(k_vec_lab[0]) * np.arccos(
        -k_vec_lab[2] / (k_magnitude * cos_chi)
    )
    theta_rad = delta_rad + tth_rad / 2
    theta = np.degrees(theta_rad)  # theta is done calculating

    sin_delta = np.sin(delta_rad)
    cos_delta = np.cos(delta_rad)
    if (np.abs(k_vec_lab[0]) < epsilon) and (np.abs(k_vec_lab[1]) < epsilon):
        phi = 0
        chi = 0
    else:
        A = k_magnitude * sin_delta
        B = -k_magnitude * cos_delta * sin_chi
        C = k_magnitude * cos_delta * sin_chi
        D = k_magnitude * sin_delta

        mat = np.array([[A, B], [C, D]])
        mat_inv = np.linalg.inv(mat)
        vector = (k_vec_lab[0], k_vec_lab[1])
        vector_rotated = mat_inv @ vector
        phi_rad = np.arctan(vector_rotated[1] / vector_rotated[0])
        phi = np.degrees(phi_rad)

    return {
        "tth": tth,
        "theta": theta,
        "phi": phi,
        "chi": chi,
        "H": H,
        "K": K,
        "L": L,
        "success": True,
        "error": None,
    }


def _calculate_angles_phi_fixed(
    k_in,
    H,
    K,
    L,
    a_star_vec_lab,
    b_star_vec_lab,
    c_star_vec_lab,
    phi,
):
    """Calculate scattering angles from HKL indices, with phi fixed."""

    # get momentum transfer vector in lab coordinates
    k_vec_lab = H * a_star_vec_lab + K * b_star_vec_lab + L * c_star_vec_lab
    k_magnitude = np.linalg.norm(k_vec_lab)
    tth_rad = 2 * np.arcsin(k_magnitude / (2 * k_in))
    tth = np.degrees(tth_rad)  # tth is done calculating

    cos_phi = np.cos(np.radians(phi))
    sin_phi = np.sin(np.radians(phi))
    vector = (k_vec_lab[0] / k_magnitude, k_vec_lab[1] / k_magnitude)
    mat = np.array([[cos_phi, -sin_phi], [sin_phi, cos_phi]])
    mat_inv = np.linalg.inv(mat)
    vector_rotated = mat_inv @ vector
    v0, v1 = vector_rotated
    v2 = -k_vec_lab[2] / k_magnitude

    chi_rad = np.arctan(v1 / v2)
    chi = np.degrees(chi_rad)

    delta = np.arcsin(v0)

    theta_rad = delta + tth_rad / 2
    theta = np.degrees(theta_rad)

    return {
        "tth": tth,
        "theta": theta,
        "phi": phi,
        "chi": chi,
        "H": H,
        "K": K,
        "L": L,
        "success": True,
        "error": None,
    }


def _calculate_angles_tth_fixed(
    k_in,
    tth,
    H=0.15,
    K=0.1,
    L=None,
    a_vec_lab=None,
    b_vec_lab=None,
    c_vec_lab=None,
    a_star_vec_lab=None,
    b_star_vec_lab=None,
    c_star_vec_lab=None,
    fixed_angle_name="chi",
    fixed_angle=0.0,
):
    """Calculate scattering angles from two of the three HKL indices, with tth fixed."""

    # initial k_vec_lab when sample has not rotated
    k_vec_lab = calculate_k_vector_in_lab(k_in, tth)

    print(f"calculated {np.dot(k_vec_lab, a_vec_lab)/np.pi/2}, H")
    H = H if H is not None else np.dot(k_vec_lab, a_vec_lab) / (2 * np.pi)
    K = K if K is not None else np.dot(k_vec_lab, b_vec_lab) / (2 * np.pi)
    L = L if L is not None else np.dot(k_vec_lab, c_vec_lab) / (2 * np.pi)

    k_vec_temp = H * a_star_vec_lab + K * b_star_vec_lab + L * c_star_vec_lab
    assert np.allclose(
        k_vec_temp, k_vec_lab
    ), f"k_vec_temp: {k_vec_temp}, k_vec_lab: {k_vec_lab}"
    calculate_angles = _calculate_angles_factory(fixed_angle_name)
    result = calculate_angles(
        k_in, H, K, L, a_star_vec_lab, b_star_vec_lab, c_star_vec_lab, fixed_angle
    )
    assert (
        np.abs(result["tth"] - tth) < 1e-6
    ), f"tth_calculated: {result['tth']}, tth_input: {tth}"
    return result


def _get_real_space_vectors(a, b, c, alpha, beta, gamma):
    """Get the real space vectors a_vec, b_vec, c_vec from the lattice parameters.
    - a_vec is by-default along x-axis (a, 0, 0)
    - b_vec is by-default (b cos gamma, b sin gamma, 0) on the x-y plane,
    - c_vec is then calculated
    The above convention defines the crystal coordinate system.

    Args:
        a, b, c (float): Lattice constants in Angstroms
        alpha, beta, gamma (float): Lattice angles in degrees

    Returns:
        a_vec, b_vec, c_vec (np.ndarray): Real space vectors
    """
    alpha_rad, beta_rad, gamma_rad = (
        np.radians(alpha),
        np.radians(beta),
        np.radians(gamma),
    )
    a_vec = np.array([a, 0, 0])
    b_vec = np.array([b * np.cos(gamma_rad), b * np.sin(gamma_rad), 0])
    c_vec_x = c * np.cos(beta_rad)
    c_vec_y = (
        c
        * (np.cos(alpha_rad) - np.cos(beta_rad) * np.cos(gamma_rad))
        / np.sin(gamma_rad)
    )
    c_vec_z = np.sqrt(c**2 - c_vec_x**2 - c_vec_y**2)
    c_vec = np.array([c_vec_x, c_vec_y, c_vec_z])
    return a_vec, b_vec, c_vec


def _get_reciprocal_space_vectors(a, b, c, alpha, beta, gamma):
    """Get the reciprocal space vectors a_star_vec, b_star_vec, c_star_vec from the lattice
    parameters, angles in degrees. These vectors are in the crystal coordinate system.
    """
    a_vec, b_vec, c_vec = _get_real_space_vectors(a, b, c, alpha, beta, gamma)
    volumn = abs(np.dot(a_vec, np.cross(b_vec, c_vec)))
    a_star_vec = 2 * np.pi * np.cross(b_vec, c_vec) / volumn
    b_star_vec = 2 * np.pi * np.cross(c_vec, a_vec) / volumn
    c_star_vec = 2 * np.pi * np.cross(a_vec, b_vec) / volumn
    return a_star_vec, b_star_vec, c_star_vec


def _get_norm_vector(h, k, l, a, b, c, alpha, beta, gamma):
    """Get the norm vector of the plane defined by the Miller indices (h, k, l)."""
    a_star_vec, b_star_vec, c_star_vec = _get_reciprocal_space_vectors(
        a, b, c, alpha, beta, gamma
    )
    norm_vec = (
        h * a_star_vec / (2 * np.pi)
        + k * b_star_vec / (2 * np.pi)
        + l * c_star_vec / (2 * np.pi)
    )
    return norm_vec


def _get_d_spacing(h, k, l, a, b, c, alpha, beta, gamma):
    """Get the d-spacing of the plane defined by the Miller indices (h, k, l)."""
    norm_vec = _get_norm_vector(h, k, l, a, b, c, alpha, beta, gamma)
    d_spacing = 1 / np.linalg.norm(norm_vec)
    return d_spacing


def _get_momentum_diffraction(h, k, l, a, b, c, alpha, beta, gamma):
    """Get the momentum transfer vector of the plane defined by the Miller indices (h, k, l)."""
    norm_vec = _get_norm_vector(h, k, l, a, b, c, alpha, beta, gamma)
    return 2 * np.pi * norm_vec


def _get_HKL_from_momentum_scattering(momentum, a_vec, b_vec, c_vec):
    """Get the HKL (r.l.u.) from the momentum transfer vector."""
    H = np.dot(momentum, a_vec) / (2 * np.pi)
    K = np.dot(momentum, b_vec) / (2 * np.pi)
    L = np.dot(momentum, c_vec) / (2 * np.pi)
    return H, K, L


def calculate_k_magnitude(k_in, tth):
    """Calculate the momentum transfer magnitude from the scattering angle."""
    return 2 * k_in * np.sin(np.radians(tth / 2.0))


def calculate_k_vector_in_lab(k_in, tth):
    """get the momentum transfer k vector in lab frame from the scattering angle tth"""
    eta = 90 - tth / 2
    eta_rad = np.radians(eta)
    k_magnitude = calculate_k_magnitude(k_in, tth)
    k_vector = k_magnitude * np.array([-np.cos(eta_rad), 0, -np.sin(eta_rad)])
    return k_vector
