"""This is a module responsible for calculating the structure factor of a given crystal
structure."""
import Dans_Diffraction as dif
import numpy as np

class StructureFactorCalculator:
    """
    A class to calculate structure factors from a CIF file.

    The typical workflow is:
    1. Instantiate the class with the CIF file path.
    2. Call initialize() to load the CIF file.
    3. Set HKL input mode and HKL values.
    4. Set energy mode and energy values.
    5. Call calculate_structure_factors() to get the results.
    """

    # Predefined X-ray sources and their wavelengths in Angstroms
    PREDEFINED_SOURCES = {
        "Cr Kα1": 2.28973,
        "Cu Kα1": 1.54059,
        "Cu Kβ": 1.39223,
        "Mo Kα1": 0.709317,
        "Ag Kα1": 0.559422,
    }

    # Default HKL list for HKL mode 1


    # Constant for converting wavelength (Angstrom) to energy (keV)
    # E(keV) = HC_E / lambda(Angstrom)
    HC_E = 12.3984198

    def __init__(self):
        """
        Initializes the StructureFactorCalculator.

        Args:
            cif_file_path (str): Path to the .cif file.
        """
        self.cif_file_path = None
        self.xtl = None
        self._initialized = False

        self._hkl_mode = None  # 0 for user-defined single HKL, 1 for list
        self._hkl_input_list = None # Stores [[h,k,l]] or list of [[h,k,l], ...]

        self._energy_mode = None  # 0 for predefined source, 1 for custom wavelength
        self._energy_kev = None

    def initialize(self, cif_file_path):
        """
        Loads the CIF file and initializes the crystal object.

        Returns:
            bool: True if initialization is successful, False otherwise.
        """
        self.cif_file_path = cif_file_path
        try:
            self.xtl = dif.Crystal(self.cif_file_path)
            self._initialized = True
            print(f"Successfully loaded CIF file: {self.cif_file_path}")
            return True
        except FileNotFoundError:
            print(f"Error: CIF file not found at '{self.cif_file_path}'. Please ensure the path is correct.")
            self._initialized = False
            return False
        except Exception as e:
            print(f"Error loading CIF file: {e}")
            self._initialized = False
            return False

    @property
    def is_initialized(self):
        return self._initialized

    @property
    def hkl_mode(self):
        return self._hkl_mode

    @hkl_mode.setter
    def hkl_mode(self, mode: int):
        if mode not in [0, 1]:
            raise ValueError("HKL input mode must be 0 (single HKL) or 1 (HKL list).")
        self._hkl_mode = mode
        self._hkl_input_list = None  # Reset HKL input when mode changes

    @property
    def hkl_input_list(self):
        return self._hkl_input_list

    @hkl_input_list.setter
    def hkl_input_list(self, hkl_list: list[list[int]]):
        if self._hkl_mode == 0:
            raise ValueError(
                "HKL input mode must be 1 to set an HKL list. Use set_hkl_input_mode(1) first."
            )
        self._hkl_input_list = hkl_list

    def _wavelength_to_kev(self, wavelength_angstrom):
        """Converts wavelength in Angstroms to energy in keV."""
        if wavelength_angstrom <= 0:
            raise ValueError("Wavelength must be positive.")
        return self.HC_E / wavelength_angstrom

    def set_hkl_input_mode(self, mode=1):
        """
        mode (int): 0 for a single user-defined HKL, 1 for a list of HKLs.
        """
        if mode not in [0, 1]:
            raise ValueError("HKL input mode must be 0 (single HKL) or 1 (HKL list).")
        self._hkl_mode = mode
        self._hkl_input_list = None # Reset HKL input when mode changes
        print(f"HKL input mode set to: {'single HKL' if mode == 0 else 'HKL list'}")

    def set_hkl_user_defined(self, h: int, k: int, l: int):
        """
        Sets a single HKL value for calculation (for HKL mode 0). h, k, l (int): Miller index h, k, l.
        """
        if self._hkl_mode != 0:
            raise TypeError("HKL input mode must be 0 to set a user-defined HKL. Use set_hkl_input_mode(0) first.")
        self._hkl_input_list = [[int(h), int(k), int(l)]]

    def set_hkl_list(self, hkl_list: list[list[int]] = None):
        """
        Sets the list of HKL values for calculation (for HKL mode 1).

        Args:
            hkl_list (list, optional): A list of HKLs. e.g. [[1,0,0], [0,1,0], [0,0,1]]

        """
        if self._hkl_mode != 1:
            raise TypeError(
                "HKL input mode must be 1 to set an HKL list. Use set_hkl_input_mode(1) first."
            )
        else:
            if not isinstance(hkl_list, list) or not all(isinstance(hkl, (list, tuple)) and len(hkl) == 3 for hkl in hkl_list):
                raise ValueError("hkl_list must be a list of lists/tuples, each containing 3 integers.")
            self._hkl_input_list = [list(map(int, hkl)) for hkl in hkl_list] # Ensure inner elements are lists of ints
            print(f"Custom HKL list set with {len(self._hkl_input_list)} entries.")

    def set_energy_mode(self, mode=0):
        """
        mode (int): 0 for predefined X-ray sources, 1 for custom wavelength.

        Raises:
            ValueError: If mode is not 0 or 1.
        """
        if mode not in [0, 1]:
            raise ValueError("Energy mode must be 0 (predefined source) or 1 (custom wavelength).")
        self._energy_mode = mode
        self._energy_kev = None # Reset energy when mode changes
        print(f"Energy input mode set to: {'predefined source' if mode == 0 else 'custom wavelength'}")

    def set_predefined_energy_source(self, source_name):
        """
        Sets the energy based on a predefined X-ray source (for energy mode 0).

        Args:
            source_name (str): The name of the predefined source (e.g., "Cu Kα1").

        Raises:
            TypeError: If energy mode is not 0.
            ValueError: If the source_name is not in PREDEFINED_SOURCES.
        """
        if self._energy_mode != 0:
            raise TypeError("Energy mode must be 0 to use a predefined source. Use set_energy_mode(0) first.")
        if source_name not in self.PREDEFINED_SOURCES:
            raise ValueError(f"Unknown predefined source: {source_name}. Available sources: {list(self.PREDEFINED_SOURCES.keys())}")

        wavelength = self.PREDEFINED_SOURCES[source_name]
        self._energy_kev = self._wavelength_to_kev(wavelength)
        print(f"Energy set from predefined source: {source_name} ({self._energy_kev:.6f} keV)")

    def set_custom_energy_wavelength(self, wavelength_angstrom):
        """
        Sets the energy based on a custom wavelength (for energy mode 1).

        Args:
            wavelength_angstrom (float): The wavelength in Angstroms.

        Raises:
            TypeError: If energy mode is not 1.
            ValueError: If wavelength is not positive.
        """
        if self._energy_mode != 1:
            raise TypeError("Energy mode must be 1 to set a custom wavelength. Use set_energy_mode(1) first.")
        if not isinstance(wavelength_angstrom, (int, float)) or wavelength_angstrom <= 0:
            raise ValueError("Wavelength must be a positive number.")

        self._energy_kev = self._wavelength_to_kev(float(wavelength_angstrom))
        print(f"Energy set from custom wavelength: {wavelength_angstrom} Å ({self._energy_kev:.6f} keV)")

    def calculate_structure_factors(self):
        """
        Calculates the structure factors based on the current settings.

        Returns:
            numpy.ndarray or None: 
                - For HKL mode 0: A NumPy array [real_part, imaginary_part].
                - For HKL mode 1: A NumPy array [real1, imag1, real2, imag2, ...].
                - Returns None if calculation cannot be performed due to missing setup or errors.
        """
        if not self.is_initialized:
            print("Error: Calculator is not initialized. Call initialize() first.")
            return None
        if self._hkl_mode is None or self._hkl_input_list is None:
            print("Error: HKL parameters are not set. Use set_hkl_input_mode() and related methods.")
            return None
        if self._energy_mode is None or self._energy_kev is None:
            print("Error: Energy parameters are not set. Use set_energy_mode() and related methods.")
            return None

        try:
            # Setup scattering: output=False suppresses console output from Dans_Diffraction
            self.xtl.Scatter.setup_scatter(scattering_type='x-ray', energy_kev=self._energy_kev, output=False)

            # Calculate structure factors
            # 'xray dispersion' typically includes anomalous scattering effects.
            F_hkl_complex_array = self.xtl.Scatter.structure_factor(
                hkl=self._hkl_input_list,
                scattering_type='xray dispersion'
            )

            if not isinstance(F_hkl_complex_array, np.ndarray):
                F_hkl_complex_array = np.array(F_hkl_complex_array)

            results = []
            if self._hkl_mode == 0: # Single HKL
                if F_hkl_complex_array.ndim == 0: # Scalar complex number for single HKL
                    F_complex = F_hkl_complex_array.item() # Get the complex number itself
                elif F_hkl_complex_array.size == 1: # Array with one complex number
                    F_complex = F_hkl_complex_array[0]
                else:
                    print(f"Error: Expected a single complex value for mode 0, but got array of shape {F_hkl_complex_array.shape}")
                    return None
                results.extend([F_complex.real, F_complex.imag])

            elif self._hkl_mode == 1: # List of HKLs
                if len(F_hkl_complex_array) == len(self._hkl_input_list):
                    for F_complex in F_hkl_complex_array:
                        results.extend([F_complex.real, F_complex.imag])
                else:
                    print(f"Error: Mismatch between number of HKL inputs ({len(self._hkl_input_list)}) and F_hkl results ({len(F_hkl_complex_array)}).")
                    return None

            return np.array(results)

        except Exception as e:
            print(f"Error during structure factor calculation: {e}")
            return None

if __name__ == '__main__':
    # This is an example of how to use the StructureFactorCalculator class.
    # You need to have a CIF file (e.g., "test.cif") in the same directory or provide the full path.
    # Create a dummy CIF file for testing if you don't have one.
    # For example, save the following as "NaCl.cif":
    """
data_NaCl
_symmetry_space_group_name_H-M   'F m -3 m'
_cell_length_a   5.6402
_cell_length_b   5.6402
_cell_length_c   5.6402
_cell_angle_alpha   90
_cell_angle_beta    90
_cell_angle_gamma   90

loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Na Na+ 0.00000 0.00000 0.00000 1.0
Cl Cl- 0.50000 0.00000 0.00000 1.0
    """
    # Ensure you have Dans_Diffraction installed: pip install Dans_Diffraction

    print("--- Example Usage of StructureFactorCalculator ---")

    # Path to your CIF file
    cif_file = "data/nacl.cif"  # Replace with your CIF file path

    # 1. Create an instance of the calculator
    calculator = StructureFactorCalculator()

    # 2. Initialize the calculator (loads CIF)
    if not calculator.initialize(cif_file_path=cif_file):
        print("Exiting due to initialization failure.")
        exit()

    # --- Scenario 1: Single HKL (1,1,1), Cu Kα1 energy ---
    print("\n--- Scenario 1: Single HKL (1,1,1), Cu Kα1 energy ---")
    calculator.set_hkl_input_mode(0)
    calculator.set_hkl_user_defined(1, 1, 1)

    calculator.set_energy_mode(0)
    calculator.set_predefined_energy_source("Cu Kα1")

    results_scenario1 = calculator.calculate_structure_factors()
    if results_scenario1 is not None:
        print(f"Results for (1,1,1) [Real(F), Imag(F)]: {results_scenario1}")
        # Example: F(111) for NaCl with Cu Kα1. Na: f=11, Cl: f=17.
        # F(111) = 4 * (f_Na - f_Cl) if origin at Na.
        # If origin is at center of symmetry, F(111) = 4 * (f_Na * exp(0) + f_Cl * exp(i*pi*(h+k+l)))
        # For hkl all odd (like 111), F = 4(f_Na - f_Cl)
        # For hkl all even, F = 4(f_Na + f_Cl)
        # For mixed parity, F = 0
        # This is a simplified view; Dans_Diffraction handles atomic form factors correctly.

    # --- Scenario 2: List of HKLs, custom energy (Mo Kα1 wavelength) ---
    print("\n--- Scenario 2: List of HKLs, custom energy (Mo Kα1 wavelength) ---")
    calculator.set_hkl_input_mode(1)
    # Using default HKL list, or provide your own:
    # calculator.set_hkl_list([[1,0,0], [2,0,0], [2,2,0]])
    calculator.set_hkl_list() # Uses default list

    calculator.set_energy_mode(1)
    mo_kalpha1_wavelength = StructureFactorCalculator.PREDEFINED_SOURCES["Mo Kα1"]
    calculator.set_custom_energy_wavelength(mo_kalpha1_wavelength)

    results_scenario2 = calculator.calculate_structure_factors()
    if results_scenario2 is not None:
        print(f"Results for HKL list (flat array [R1,I1, R2,I2 ...]):")
        # Reshape for easier viewing if needed
        num_hkls = len(calculator._hkl_input_list)
        reshaped_results = results_scenario2.reshape((num_hkls, 2))
        for i, hkl in enumerate(calculator._hkl_input_list):
            print(f"  HKL {tuple(hkl)}: F_real = {reshaped_results[i,0]:.4f}, F_imag = {reshaped_results[i,1]:.4f}, |F|^2 = {reshaped_results[i,0]**2 + reshaped_results[i,1]**2:.4f}")

    # --- Scenario 3: Invalid usage (e.g., source not found) ---
    print("\n--- Scenario 3: Error Handling Example ---")
    calculator.set_energy_mode(0)
    try:
        calculator.set_predefined_energy_source("NonExistentSource")
    except ValueError as e:
        print(f"Caught expected error: {e}")

    # --- Scenario 4: Using a different HKL from default list for single HKL mode ---
    print("\n--- Scenario 4: Single HKL (2,0,0), Cu Kα1 energy ---")
    calculator.set_hkl_input_mode(0)
    calculator.set_hkl_user_defined(2,0,0) # (2,0,0) is an all-even reflection for NaCl

    calculator.set_energy_mode(0) # Already set to Cu Kα1 from scenario 1 if energy settings were not changed
    calculator.set_predefined_energy_source("Cu Kα1") # Re-set to be sure

    results_scenario4 = calculator.calculate_structure_factors()
    if results_scenario4 is not None:
        print(f"Results for (2,0,0) [Real(F), Imag(F)]: {results_scenario4}")
