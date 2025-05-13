#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# import partial
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from functools import partial
import numpy as np
from scipy.optimize import fsolve
from packages.classes import Lab
from packages.utils import angle_to_matrix


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
        print(f"k_in: {self.k_in}, k_magnitude: {k_magnitude}")
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


def calculate_tth_from_k_magnitude(k_in, k_magnitude):
    """calculate the scattering angle tth from the momentum transfer magnitude"""
    return 2 * np.degrees(np.arcsin(k_magnitude / (2 * k_in)))


def calculate_k_vector_in_lab(k_in, tth):
    """get the momentum transfer k vector in lab frame from the scattering angle tth"""
    eta = 90 - tth / 2
    eta_rad = np.radians(eta)
    k_magnitude = calculate_k_magnitude(k_in, tth)
    k_vector = k_magnitude * np.array([-np.cos(eta_rad), 0, -np.sin(eta_rad)])
    return k_vector


def derivative(fun, x, delta_x=1e-6):
    """calculate the derivative of the function fun at the point x"""
    return (fun(x + delta_x) - fun(x - delta_x)) / (2 * delta_x)


def process_angle(angle):
    """process the angle to be in the range of (-180, 180]"""
    angle = angle % 360
    if angle > 180:
        angle -= 360
    return angle


def _calculate_angles_factory(fixed_angle_name):
    if fixed_angle_name == "chi":
        return _calculate_angles_chi_fixed
    elif fixed_angle_name == "phi":
        return _calculate_angles_phi_fixed
    elif fixed_angle_name == "tth":
        return _calculate_angles_tth_fixed


def _calculate_angles_tth_fixed(
    k_in,
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
    H=0.15,
    K=0.1,
    L=None,
    fixed_angle_name="chi",
    fixed_angle=0.0,
    target_objective=1e-4,
    num_steps=1000,
    number_batch=10,
    learning_rate=100,
):
    """Calculate scattering angles from two of the three HKL indices, with tth fixed."""

    # initial k_vec_lab when sample has not rotated
    k_magnitude_target = calculate_k_magnitude(k_in, tth)
    lab = Lab()
    lab.initialize(a, b, c, alpha, beta, gamma, roll, pitch, yaw, 0, 0, 0)
    a_star_vec_lab, b_star_vec_lab, c_star_vec_lab = lab.get_reciprocal_space_vectors()

    # Define which index is None and will be solved for
    index_to_solve = None
    if H is None:
        index_to_solve = "H"
    elif K is None:
        index_to_solve = "K"
    elif L is None:
        index_to_solve = "L"

    def fun_to_solve(momentum):
        h_val = momentum if index_to_solve == "H" else H
        k_val = momentum if index_to_solve == "K" else K
        l_val = momentum if index_to_solve == "L" else L
        k = h_val * a_star_vec_lab + k_val * b_star_vec_lab + l_val * c_star_vec_lab
        k_magnitude = np.linalg.norm(k)
        return k_magnitude - k_magnitude_target

    momentum = fsolve(fun_to_solve, -0.4)

    # Update the appropriate index
    if index_to_solve == "H":
        H = momentum[0]
    elif index_to_solve == "K":
        K = momentum[0]
    elif index_to_solve == "L":
        L = momentum[0]

    print(f"H: {H}, K: {K}, L: {L}")
    calculate_angles = _calculate_angles_factory(fixed_angle_name)

    result = calculate_angles(
        k_in,
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
    return result


def _calculate_angles_chi_fixed(
    k_in,
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
    chi_fixed,
    target_objective=1e-4,
    num_steps=1000,
    number_batch=10,
    learning_rate=100,
):
    """gradient decent trial function"""

    def objective_function(k_cal, k_target):
        """objective function for gradient decent"""
        return np.linalg.norm(k_cal - k_target)

    def get_k_cal(lab, theta_, phi_, chi_):
        lab.rotate(theta_, phi_, chi_)
        a_star_vec, b_star_vec, c_star_vec = lab.get_reciprocal_space_vectors()
        k_cal = H * a_star_vec + K * b_star_vec + L * c_star_vec
        return k_cal

    theta_best_list = [0] * number_batch
    phi_best_list = [0] * number_batch
    for batch in range(number_batch):
        lab = Lab()
        theta = np.random.uniform(0, 180)
        phi = np.random.uniform(0, 180)

        lab.initialize(
            a, b, c, alpha, beta, gamma, roll, pitch, yaw, theta, phi, chi_fixed
        )

        k_cal = get_k_cal(lab, theta, phi, chi_fixed)
        k_magnitude = np.linalg.norm(k_cal)
        tth = calculate_tth_from_k_magnitude(k_in, k_magnitude)
        print(f"tth: {tth}")
        k_target = calculate_k_vector_in_lab(k_in, tth)
        objective = objective_function(k_cal, k_target)
        for i in range(num_steps):
            step_size = objective * learning_rate
            theta_new = theta + np.random.uniform(-step_size, step_size)
            phi_new = phi + np.random.uniform(-step_size, step_size)
            k_cal = get_k_cal(lab, theta_new, phi_new, chi_fixed)
            objective_new = objective_function(k_cal, k_target)
            if objective_new < objective:
                theta = theta_new
                phi = phi_new
                objective = objective_new
            if objective < target_objective:
                break
        # Normalize angles to (0, 360) range
        theta = process_angle(theta)
        phi = process_angle(phi)

        theta_best_list[batch] = theta
        phi_best_list[batch] = phi

    # round up to 0.1, discard duplicates, phi and theta should match the order of the list
    theta_best_list = np.round(theta_best_list, 1)
    phi_best_list = np.round(phi_best_list, 1)
    theta_result = []
    phi_result = []
    tth_result = []
    chi_result = []
    for theta, phi in zip(theta_best_list, phi_best_list):
        if theta not in theta_result:
            theta_result.append(theta)
            phi_result.append(phi)
            tth_result.append(process_angle(tth))
            chi_result.append(chi_fixed)

    return tth_result, theta_result, phi_result, chi_result


def _calculate_angles_phi_fixed(
    k_in,
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
    phi_fixed,
    target_objective=1e-4,
    num_steps=1000,
    number_batch=10,
    learning_rate=100,
):
    """gradient decent trial function with phi fixed"""

    def objective_function(k_cal, k_target):
        """objective function for gradient decent"""
        return np.linalg.norm(k_cal - k_target)

    def get_k_cal(lab, theta_, phi_, chi_):
        lab.rotate(theta_, phi_, chi_)
        a_star_vec, b_star_vec, c_star_vec = lab.get_reciprocal_space_vectors()
        k_cal = H * a_star_vec + K * b_star_vec + L * c_star_vec
        return k_cal

    theta_best_list = [0] * number_batch
    chi_best_list = [0] * number_batch
    for batch in range(number_batch):
        lab = Lab()
        theta = np.random.uniform(0, 180)
        chi = np.random.uniform(0, 180)

        lab.initialize(
            a, b, c, alpha, beta, gamma, roll, pitch, yaw, theta, phi_fixed, chi
        )

        k_cal = get_k_cal(lab, theta, phi_fixed, chi)
        k_magnitude = np.linalg.norm(k_cal)
        tth = calculate_tth_from_k_magnitude(k_in, k_magnitude)
        k_target = calculate_k_vector_in_lab(k_in, tth)
        objective = objective_function(k_cal, k_target)
        for i in range(num_steps):
            step_size = objective * learning_rate
            theta_new = theta + np.random.uniform(-step_size, step_size)
            chi_new = chi + np.random.uniform(-step_size, step_size)
            k_cal = get_k_cal(lab, theta_new, phi_fixed, chi_new)
            objective_new = objective_function(k_cal, k_target)
            if objective_new < objective:
                theta = theta_new
                chi = chi_new
                objective = objective_new
            if objective < target_objective:
                break
        # Normalize angles to (0, 360) range
        theta = process_angle(theta)
        chi = process_angle(chi)
        theta_best_list[batch] = theta
        chi_best_list[batch] = chi

    # round up to 0.1, discard duplicates, theta and chi should match the order of the list
    theta_best_list = np.round(theta_best_list, 1)
    chi_best_list = np.round(chi_best_list, 1)
    theta_result = []
    chi_result = []
    tth_result = []
    phi_result = []
    for theta, chi in zip(theta_best_list, chi_best_list):
        if theta not in theta_result:
            theta_result.append(theta)
            chi_result.append(chi)
            tth_result.append(process_angle(tth))
            phi_result.append(phi_fixed)

    return tth_result, theta_result, phi_result, chi_result


if __name__ == "__main__":
    params = {
        "k_in": 0.48143441803485754,
        "tth": 131.98,
        "a": 1,
        "b": 2,
        "c": 3,
        "alpha": 90,
        "beta": 90,
        "gamma": 90,
        "roll": 0,
        "pitch": 0,
        "yaw": 0,
        "fixed_angle_name": "chi",
        "fixed_angle": 0,
        "H": -0.0569,
        "K": 0,
        "L": None,  # L = -0.3837
        "target_objective": 1e-4,
        "num_steps": 1000,
        "number_batch": 5,
    }
    tth_results, theta_results, phi_results, chi_results = _calculate_angles_tth_fixed(
        **params
    )
    print(f"FINAL theta: {theta_results}, phi: {phi_results}, chi: {chi_results}")
    print(f"tth: {tth_results}")
