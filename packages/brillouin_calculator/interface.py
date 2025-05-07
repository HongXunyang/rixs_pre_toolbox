#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# import partial
from functools import partial
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
        self.a_vec, self.b_vec, self.c_vec = _get_real_space_vectors(
            self.a, self.b, self.c, self.alpha, self.beta, self.gamma
        )
        self.lattice_type = "cubic"
        self.a_star_vec, self.b_star_vec, self.c_star_vec = (
            _get_reciprocal_space_vectors(
                self.a, self.b, self.c, self.alpha, self.beta, self.gamma
            )
        )
        self.e_H, self.e_K, self.e_L = (
            np.array([1, 0, 0]),
            np.array([0, 1, 0]),
            np.array([0, 0, 1]),
        )

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
            lattice_type (str): Lattice type, one of "cubic"... (MARKED: TO BE EXTENDED)
        Returns:
            bool: True if initialization was successful
        """
        if params["e_H"] is None or params["e_K"] is None or params["e_L"] is None:
            params["e_H"] = np.array([1, 0, 0])
            params["e_K"] = np.array([0, 1, 0])
            params["e_L"] = np.array([0, 0, 1])
        try:
            # Store parameters
            self.a = params["a"]
            self.b = params["b"]
            self.c = params["c"]
            self.alpha = params["alpha"]
            self.beta = params["beta"]
            self.gamma = params["gamma"]
            self.energy = params["energy"]
            self.e_H, self.e_K, self.e_L = params["e_H"], params["e_K"], params["e_L"]

            # Calculate wavelength and wavevector
            self.lambda_A = (
                (self.hPlanck * self.c_light) / (self.energy * self.e) * 1e10
            )
            self.k_in = 1.0 / self.lambda_A

            # Calculate reciprocal lattice
            self._calculate_reciprocal_lattice()

            # Initialize visualizer
            # self.visualizer = BrillouinVisualizer(self.reciprocal_lattice)

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
        self.a_star_vec, self.b_star_vec, self.c_star_vec = (
            _get_reciprocal_space_vectors(
                self.a, self.b, self.c, self.alpha, self.beta, self.gamma
            )
        )
        self.a_vec, self.b_vec, self.c_vec = _get_real_space_vectors(
            self.a, self.b, self.c, self.alpha, self.beta, self.gamma
        )

    def _Q_magnitude(self, tth):
        return 2.0 * self.k_in * np.sin(np.radians(tth / 2.0))

    def calculate_hkl(self, tth, theta, phi, chi):
        """Calculate HKL indices from scattering angles.

        Args:
            tth (float): Scattering angle in degrees
            theta (float): Sample theta rotation in degrees
            phi (float): Sample phi rotation in degrees
            chi (float): Sample chi rotation in degrees

        Returns:
            dict: Dictionary containing h, k, l (crystal coordinates) indices and Q magnitude
        """
        if not self.is_initialized():
            raise ValueError("Calculator not initialized")

        # Calculate Q magnitude
        Q = self._Q_magnitude(tth)
        # Calculate delta = theta - (tth/2)
        delta = theta - (tth / 2.0)

        # Build rotation matrices for chi and phi
        chi_rad = np.radians(chi)
        phi_rad = np.radians(phi)

        # the rotation of the sample is defined to be counter-clockwise
        chi_mat_sample = np.array(
            [
                [1, 0, 0],
                [0, np.cos(chi_rad), -np.sin(chi_rad)],
                [0, np.sin(chi_rad), np.cos(chi_rad)],
            ]
        )
        # the rotation of the beam relative to the sample is therefore clockwise
        chi_mat_beam = chi_mat_sample.T

        # the rotation of the sample is defined to be counter-clockwise
        phi_mat_sample = np.array(
            [
                [np.cos(phi_rad), -np.sin(phi_rad), 0],
                [np.sin(phi_rad), np.cos(phi_rad), 0],
                [0, 0, 1],
            ]
        )
        # the rotation of the beam relative to the sample is therefore clockwise
        phi_mat_beam = phi_mat_sample.T

        # Q vector in "theta-2theta coordinates"
        sin_delta = np.sin(np.radians(delta))
        cos_delta = np.cos(np.radians(delta))
        Q_th2th = np.array([Q * sin_delta, 0.0, -Q * cos_delta])

        # Transform to sample coordinates
        q_sample = phi_mat_beam @ chi_mat_beam @ Q_th2th
        H_lab = q_sample[0] * self.a
        K_lab = q_sample[1] * self.b
        L_lab = q_sample[2] * self.c
        H_crystal, K_crystal, L_crystal = _lab_to_crystal_coordinate(
            self.a,
            self.b,
            self.c,
            H_lab,
            K_lab,
            L_lab,
            self.e_H,
            self.e_K,
            self.e_L,
        )
        return {
            "H": H_crystal,
            "K": K_crystal,
            "L": L_crystal,
            "tth": tth,
            "theta": theta,
            "phi": phi,
            "chi": chi,
            "success": True,
            "error": None,
        }

    def calculate_angles(
        self,
        H_crystal,
        K_crystal,
        L_crystal,
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
        return calculate_angles(
            self.a,
            self.b,
            self.c,
            self.k_in,
            H_crystal,
            K_crystal,
            L_crystal,
            self.e_H,
            self.e_K,
            self.e_L,
            fixed_angle,
        )

    def calculate_angles_tth_fixed(
        self,
        tth,
        H_crystal=0.15,
        K_crystal=0.1,
        L_crystal=None,
        fixed_angle_name="chi",
        fixed_angle=0.0,
    ):
        """Calculate scattering angles from two of the three HKL indices, with tth fixed."""
        if not self.is_initialized():
            raise ValueError("Calculator not initialized")

        e_H, e_K, e_L = self.e_H, self.e_K, self.e_L
        a, b, c = self.a, self.b, self.c

        H_crystal_temp = H_crystal if H_crystal is not None else 0.0
        K_crystal_temp = K_crystal if K_crystal is not None else 0.0
        L_crystal_temp = L_crystal if L_crystal is not None else 0.0

        kh_crystal = H_crystal_temp / a
        kk_crystal = K_crystal_temp / b
        kl_crystal = L_crystal_temp / c

        Q_magnitude = self._Q_magnitude(tth)
        remainder = -np.sqrt(
            Q_magnitude**2 - kh_crystal**2 - kk_crystal**2 - kl_crystal**2
        )

        kh_crystal = kh_crystal if H_crystal is not None else remainder
        kk_crystal = kk_crystal if K_crystal is not None else remainder
        kl_crystal = kl_crystal if L_crystal is not None else remainder

        H_crystal, K_crystal, L_crystal = kh_crystal * a, kk_crystal * b, kl_crystal * c
        H_lab, K_lab, L_lab = _crystal_to_lab_coordinate(
            a, b, c, H_crystal, K_crystal, L_crystal, e_H, e_K, e_L
        )
        calculate_angles = _calculate_angles_factory(fixed_angle_name)
        result = calculate_angles(
            a, b, c, self.k_in, H_lab, K_lab, L_lab, e_H, e_K, e_L, fixed_angle
        )
        assert np.abs(result["tth"] - tth) < 1e-6
        print(f"input tth: {tth}, output tth: {result['tth']}")
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
        return {
            "a": self.a,
            "b": self.b,
            "c": self.c,
            "alpha": self.alpha,
            "beta": self.beta,
            "gamma": self.gamma,
        }


def _calculate_angles_factory(fixed_angle_name):
    if fixed_angle_name == "chi":
        return _calculate_angles_chi_fixed
    elif fixed_angle_name == "phi":
        return _calculate_angles_phi_fixed


def _calculate_angles_chi_fixed(
    a, b, c, k_in, H_crystal, K_crystal, L_crystal, e_H, e_K, e_L, chi
):
    """Calculate scattering angles from HKL (crystal coordinates) indices, with chi fixed."""

    # convert crystal coordinates to lab coordinates
    epsilon = 1e-10
    H_lab, K_lab, L_lab = _crystal_to_lab_coordinate(
        a, b, c, H_crystal, K_crystal, L_crystal, e_H, e_K, e_L
    )

    kh_lab = H_lab / a  # k_h has a unit of 1/Å
    kk_lab = K_lab / b  # k_k has a unit of 1/Å
    kl_lab = L_lab / c  # k_l has a unit of 1/Å
    cos_chi = np.cos(np.radians(chi))
    sin_chi = np.sin(np.radians(chi))
    Q_magnitude = np.sqrt(kh_lab**2 + kk_lab**2 + kl_lab**2)
    tth_rad = 2 * np.arcsin(Q_magnitude / (2 * k_in))
    tth = np.degrees(tth_rad)

    delta_rad = np.arccos(-kl_lab / (Q_magnitude * cos_chi))

    print(
        f"-kl_lab: {kl_lab}, Q_magnitude: {Q_magnitude}, cos_chi: {cos_chi}, -kl/(Q*cos_chi): {-kl_lab/(Q_magnitude*cos_chi)}"
    )
    theta_rad = delta_rad + tth_rad / 2
    theta = np.degrees(theta_rad)
    sin_delta = np.sin(delta_rad)
    cos_delta = np.cos(delta_rad)
    if (np.abs(H_lab) < epsilon) and (np.abs(K_lab) < epsilon):
        phi = 0
        chi = 0
    else:
        A = Q_magnitude * sin_delta
        B = -Q_magnitude * cos_delta * sin_chi
        C = Q_magnitude * cos_delta * sin_chi
        D = Q_magnitude * sin_delta

        mat = np.array([[A, B], [C, D]])
        print(f"mat: {mat}")
        mat_inv = np.linalg.inv(mat)
        vector = (kh_lab, kk_lab)
        vector_rotated = mat_inv @ vector
        phi_rad = np.arctan(vector_rotated[1] / vector_rotated[0])
        phi = np.degrees(phi_rad)

    return {
        "tth": tth,
        "theta": theta,
        "phi": phi,
        "chi": chi,
        "H": H_crystal,
        "K": K_crystal,
        "L": L_crystal,
        "success": True,
        "error": None,
    }


def _calculate_angles_phi_fixed(
    a, b, c, k_in, H_crystal, K_crystal, L_crystal, e_H, e_K, e_L, phi
):
    """Calculate scattering angles from HKL indices, with phi fixed."""

    # convert crystal coordinates to lab coordinates
    h_lab, k_lab, l_lab = _crystal_to_lab_coordinate(
        a, b, c, H_crystal, K_crystal, L_crystal, e_H, e_K, e_L
    )
    kh_lab = h_lab / a  # k_h has a unit of 1/Å
    kk_lab = k_lab / b  # k_k has a unit of 1/Å
    kl_lab = l_lab / c  # k_l has a unit of 1/Å
    cos_phi = np.cos(np.radians(phi))
    sin_phi = np.sin(np.radians(phi))

    Q_magnitude = np.sqrt(kh_lab**2 + kk_lab**2 + kl_lab**2)
    tth_rad = 2 * np.arcsin(Q_magnitude / (2 * k_in))
    tth = np.degrees(tth_rad)

    vector = (kh_lab / Q_magnitude, kk_lab / Q_magnitude)
    mat = np.array([[cos_phi, -sin_phi], [sin_phi, cos_phi]])
    mat_inv = np.linalg.inv(mat)
    vector_rotated = mat_inv @ vector
    v0, v1 = vector_rotated
    v2 = -kl_lab / Q_magnitude

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
        "H": H_crystal,
        "K": K_crystal,
        "L": L_crystal,
        "success": True,
        "error": None,
    }


def _lab_to_crystal_coordinate(a, b, c, H_lab, K_lab, L_lab, e_H, e_K, e_L):
    """Convert the momentum transfer from lab coordinates to crystal coordinates, in the unit of
    r.l.u. Convention: variable name starting with `k` is in the unit of 1/Å, and variable name of
    `h`, `k`, `l` is in the unit of r.l.u.

    Args:
        h_lab, k_lab, l_lab (float): Momentum transfer in lab coordinates, unit: r.l.u.
        e_H, e_K, e_L (np.ndarray): Unit vectors of the crystal coordinates in the lab coordinates

    Returns:
        h_crystal, k_crystal, l_crystal (float): Momentum transfer in crystal coordinates, unit: r.l.u.

    """
    # convert momentum transfer to 1/Å
    kh_lab = H_lab / a
    kk_lab = K_lab / b
    kl_lab = L_lab / c

    rotation_matrix = np.array([e_H, e_K, e_L]).T
    rotation_matrix_inv = np.linalg.inv(rotation_matrix)

    k_lab = np.array([kh_lab, kk_lab, kl_lab])
    k_crystal = rotation_matrix_inv @ k_lab
    H_crystal = k_crystal[0] * a  # convert back to r.l.u.
    K_crystal = k_crystal[1] * b  # unit: r.l.u.
    L_crystal = k_crystal[2] * c  # unit: r.l.u.
    return H_crystal, K_crystal, L_crystal


def _crystal_to_lab_coordinate(a, b, c, H_crystal, K_crystal, L_crystal, e_H, e_K, e_L):
    """Convert the momentum transfer from crystal coordinates to lab coordinates, in the unit of
    r.l.u. Convention: variable name starting with `k` is in the unit of 1/Å, and variable name of
    `h`, `k`, `l` is in the unit of r.l.u.
    """
    kh_crystal = H_crystal / a
    kk_crystal = K_crystal / b
    kl_crystal = L_crystal / c

    rotation_matrix = np.array([e_H, e_K, e_L]).T
    k_crystal = np.array([kh_crystal, kk_crystal, kl_crystal])
    k_lab = rotation_matrix @ k_crystal

    H_lab = k_lab[0] * a
    K_lab = k_lab[1] * b
    L_lab = k_lab[2] * c
    return H_lab, K_lab, L_lab


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
    volumn = np.dot(a_vec, np.cross(b_vec, c_vec))
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


def _get_momentum_transfer(h, k, l, a, b, c, alpha, beta, gamma):
    """Get the momentum transfer vector of the plane defined by the Miller indices (h, k, l)."""
    norm_vec = _get_norm_vector(h, k, l, a, b, c, alpha, beta, gamma)
    return 2 * np.pi * norm_vec
