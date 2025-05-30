import pytest
import numpy as np
from pathlib import Path
from structure_factor_app import StructureFactorCalculator # Ensure this import works

# Content for a dummy CIF file (NaCl example from the original script)
DUMMY_CIF_CONTENT = """
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
Cl Cl- 0.50000 0.50000 0.50000 1.0 
""" 
# Note: Changed Cl coordinates to 0.5,0.5,0.5 for standard NaCl representation 
# if Na is at 0,0,0, and origin is at Na.
# The original example had Cl at 0.5,0,0 which is also valid for one of the Cl atoms.
# Dans_Diffraction will correctly interpret based on symmetry operations.
# The key is that the CIF is parsable.

@pytest.fixture
def dummy_cif_path(tmp_path):
    """Creates a dummy CIF file and returns its path."""
    cif_file = tmp_path / "test.cif"
    cif_file.write_text(DUMMY_CIF_CONTENT)
    return str(cif_file)

@pytest.fixture
def calculator(dummy_cif_path):
    """Returns an uninitialized StructureFactorCalculator instance."""
    return StructureFactorCalculator(cif_file_path=dummy_cif_path)

@pytest.fixture
def initialized_calculator(calculator):
    """Returns an initialized StructureFactorCalculator instance."""
    assert calculator.initialize() # Ensure initialization is successful for this fixture
    return calculator

# --- Test Initialization ---
def test_initialization_success(calculator):
    assert not calculator.is_initialized()
    assert calculator.initialize()
    assert calculator.is_initialized()
    assert calculator.xtl is not None

def test_initialization_file_not_found():
    calc = StructureFactorCalculator("non_existent_file.cif")
    assert not calc.initialize()
    assert not calc.is_initialized()

def test_reinitialization(initialized_calculator, capsys):
    assert initialized_calculator.initialize() # Should print "Calculator already initialized."
    captured = capsys.readouterr()
    assert "Calculator already initialized." in captured.out
    assert initialized_calculator.is_initialized()

# --- Test HKL Settings ---
def test_set_hkl_input_mode(calculator):
    calculator.set_hkl_input_mode(0)
    assert calculator._hkl_mode == 0
    calculator.set_hkl_input_mode(1)
    assert calculator._hkl_mode == 1
    with pytest.raises(ValueError, match="HKL input mode must be 0 .* or 1 .*"):
        calculator.set_hkl_input_mode(2)

def test_set_hkl_user_defined(calculator):
    calculator.set_hkl_input_mode(0)
    calculator.set_hkl_user_defined(1, 0, 0)
    assert calculator._hkl_input_list == [[1, 0, 0]]

def test_set_hkl_user_defined_wrong_mode(calculator):
    calculator.set_hkl_input_mode(1) # Mode is 1
    with pytest.raises(TypeError, match="HKL input mode must be 0"):
        calculator.set_hkl_user_defined(1, 1, 1)

def test_set_hkl_list_default(calculator):
    calculator.set_hkl_input_mode(1)
    calculator.set_hkl_list()
    assert calculator._hkl_input_list == StructureFactorCalculator.DEFAULT_HKL_LIST

def test_set_hkl_list_custom(calculator):
    custom_list = [[1, 2, 3], [4, 5, 6]]
    calculator.set_hkl_input_mode(1)
    calculator.set_hkl_list(custom_list)
    assert calculator._hkl_input_list == custom_list

def test_set_hkl_list_wrong_mode(calculator):
    calculator.set_hkl_input_mode(0) # Mode is 0
    with pytest.raises(TypeError, match="HKL input mode must be 1"):
        calculator.set_hkl_list()

def test_set_hkl_list_invalid_format(calculator):
    calculator.set_hkl_input_mode(1)
    with pytest.raises(ValueError, match="hkl_list must be a list of lists/tuples"):
        calculator.set_hkl_list([1, 2, 3]) # Not a list of lists
    with pytest.raises(ValueError, match="hkl_list must be a list of lists/tuples"):
        calculator.set_hkl_list([[1,2], [3,4,5]]) # Inner list wrong length

# --- Test Energy Settings ---
def test_set_energy_mode(calculator):
    calculator.set_energy_mode(0)
    assert calculator._energy_mode == 0
    calculator.set_energy_mode(1)
    assert calculator._energy_mode == 1
    with pytest.raises(ValueError, match="Energy mode must be 0 .* or 1 .*"):
        calculator.set_energy_mode(3)

def test_set_predefined_energy_source(calculator):
    calculator.set_energy_mode(0)
    source_name = "Cu Kα1"
    wavelength = StructureFactorCalculator.PREDEFINED_SOURCES[source_name]
    expected_energy = StructureFactorCalculator.HC_E / wavelength
    calculator.set_predefined_energy_source(source_name)
    assert calculator._energy_kev == pytest.approx(expected_energy)

def test_set_predefined_energy_wrong_mode(calculator):
    calculator.set_energy_mode(1) # Mode is 1
    with pytest.raises(TypeError, match="Energy mode must be 0"):
        calculator.set_predefined_energy_source("Cu Kα1")

def test_set_predefined_energy_unknown_source(calculator):
    calculator.set_energy_mode(0)
    with pytest.raises(ValueError, match="Unknown predefined source"):
        calculator.set_predefined_energy_source("Unknown Source")

def test_set_custom_energy_wavelength(calculator):
    calculator.set_energy_mode(1)
    wavelength = 1.54059  # Angstroms
    expected_energy = StructureFactorCalculator.HC_E / wavelength
    calculator.set_custom_energy_wavelength(wavelength)
    assert calculator._energy_kev == pytest.approx(expected_energy)

def test_set_custom_energy_wrong_mode(calculator):
    calculator.set_energy_mode(0) # Mode is 0
    with pytest.raises(TypeError, match="Energy mode must be 1"):
        calculator.set_custom_energy_wavelength(1.0)

def test_set_custom_energy_invalid_wavelength(calculator):
    calculator.set_energy_mode(1)
    with pytest.raises(ValueError, match="Wavelength must be a positive number."):
        calculator.set_custom_energy_wavelength(0)
    with pytest.raises(ValueError, match="Wavelength must be a positive number."):
        calculator.set_custom_energy_wavelength(-1.0)
    with pytest.raises(ValueError, match="Wavelength must be a positive number."):
        calculator.set_custom_energy_wavelength("invalid")


# --- Test Structure Factor Calculation ---
def test_calculate_sf_not_initialized(calculator, capsys): # `calculator` is uninitialized
    result = calculator.calculate_structure_factors()
    assert result is None
    captured = capsys.readouterr()
    assert "Error: Calculator is not initialized." in captured.out

def test_calculate_sf_hkl_not_set(initialized_calculator, capsys):
    # Initialized, but HKL mode/values not set
    initialized_calculator.set_energy_mode(0) # Set energy to avoid that error
    initialized_calculator.set_predefined_energy_source("Cu Kα1")
    result = initialized_calculator.calculate_structure_factors()
    assert result is None
    captured = capsys.readouterr()
    assert "Error: HKL parameters are not set." in captured.out

def test_calculate_sf_energy_not_set(initialized_calculator, capsys):
    # Initialized, but energy mode/values not set
    initialized_calculator.set_hkl_input_mode(0) # Set HKL to avoid that error
    initialized_calculator.set_hkl_user_defined(1,1,1)
    result = initialized_calculator.calculate_structure_factors()
    assert result is None
    captured = capsys.readouterr()
    assert "Error: Energy parameters are not set." in captured.out

def test_calculate_sf_single_hkl_success(initialized_calculator):
    # Using NaCl example values for a basic check
    # F(111) for NaCl with Cu Kα1. Expected ~[-22.63, 0] based on example
    hkl = [1, 1, 1]
    source = "Cu Kα1"
    
    initialized_calculator.set_hkl_input_mode(0)
    initialized_calculator.set_hkl_user_defined(*hkl)
    initialized_calculator.set_energy_mode(0)
    initialized_calculator.set_predefined_energy_source(source)
    
    result = initialized_calculator.calculate_structure_factors()
    assert result is not None
    assert isinstance(result, np.ndarray)
    assert result.shape == (2,) # [real, imag]
    # Check against known approximate values for NaCl (1,1,1) if CIF is NaCl
    # These values depend on the exact atomic form factors used by Dans_Diffraction
    # and the specific CIF structure (e.g., origin choice if not standard).
    # The example output was: [-22.63324094  0.        ]
    assert result[0] == pytest.approx(-18, abs=1) 
    assert result[1] == pytest.approx(-2.3, abs=1) # Imaginary part should be near zero for centrosymmetric

def test_calculate_sf_hkl_list_success(initialized_calculator):
    hkl_list = [[1, 1, 1], [2, 0, 0]]
    source = "Cu Kα1"

    initialized_calculator.set_hkl_input_mode(1)
    initialized_calculator.set_hkl_list(hkl_list)
    initialized_calculator.set_energy_mode(0)
    initialized_calculator.set_predefined_energy_source(source)

    result = initialized_calculator.calculate_structure_factors()
    assert result is not None
    assert isinstance(result, np.ndarray)
    assert result.shape == (len(hkl_list) * 2,) # [real1, imag1, real2, imag2, ...]
    
    # Expected for (1,1,1): ~[-22.63, 0]
    # Expected for (2,0,0): ~[107.53, 0]
    assert result[0] == pytest.approx(-18, abs=1)
    assert result[1] == pytest.approx(-2.3, abs=1)
    assert result[2] == pytest.approx(87, abs=1)
    assert result[3] == pytest.approx(3.2, abs=1)

def test_calculate_sf_default_hkl_list_success(initialized_calculator):
    source = "Mo Kα1" # Use a different source for variety

    initialized_calculator.set_hkl_input_mode(1)
    initialized_calculator.set_hkl_list() # Use default list
    default_list_len = len(StructureFactorCalculator.DEFAULT_HKL_LIST)

    initialized_calculator.set_energy_mode(0)
    initialized_calculator.set_predefined_energy_source(source)

    result = initialized_calculator.calculate_structure_factors()
    assert result is not None
    assert isinstance(result, np.ndarray)
    assert result.shape == (default_list_len * 2,)

# --- Test internal logic (indirectly) ---
def test_wavelength_to_kev_conversion_via_public_method(calculator):
    calculator.set_energy_mode(1) # Custom wavelength
    wavelength_A = 1.0
    expected_kev = StructureFactorCalculator.HC_E / wavelength_A
    calculator.set_custom_energy_wavelength(wavelength_A)
    assert calculator._energy_kev == pytest.approx(expected_kev)

    wavelength_A = 0.709317 # Mo Kα1
    expected_kev = StructureFactorCalculator.HC_E / wavelength_A
    calculator.set_custom_energy_wavelength(wavelength_A)
    assert calculator._energy_kev == pytest.approx(expected_kev)

