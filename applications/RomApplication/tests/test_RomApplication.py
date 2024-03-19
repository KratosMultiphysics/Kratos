# import Kratos
import KratosMultiphysics

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
from test_fluid_rom import TestFluidRom
from test_thermal_rom import TestThermalRom
from test_structural_rom import TestStructuralRom
from test_randomized_singular_value_decomposition import TestRandomizedSVD
from test_empirical_cubature_method import TestEmpiricalCubatureMethod
from test_calculate_rom_basis_output_process_json import TestCalculateRomBasisOutputProcessJSON
from test_calculate_rom_basis_output_process_numpy import TestCalculateRomBasisOutputProcessNumpy
from test_compressible_potiential_rom import TestCompressiblePotentialRom
from test_fluid_lspg_rom import TestFluidLSPGRom
from test_thermal_lspg_rom import TestThermalLSPGRom
from test_structural_lspg_rom import TestStructuralLSPGRom
from test_fluid_pg_rom import TestFluidPGRom
from test_thermal_pg_rom import TestThermalPGRom
from test_structural_pg_rom import TestStructuralPGRom
from test_monotonicity_preserving_rom import TestMonotonicityPreservingRom

def AssembleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should pupulate the suites:
    "small", "nightly" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''

    suites = KratosUnittest.KratosSuites

    # Create a test suit with the selected tests (Small tests):
    # smallSuite will contain the following tests:
    # - testSmallExample
    smallSuite = suites['small']
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestFluidRom]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestThermalRom]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestStructuralRom]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCalculateRomBasisOutputProcessJSON]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCalculateRomBasisOutputProcessNumpy]))
    smallSuite.addTest(TestRandomizedSVD('test_radomized_svd'))
    smallSuite.addTest(TestEmpiricalCubatureMethod('test_empirical_cubature_method'))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCompressiblePotentialRom]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestFluidLSPGRom]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestThermalLSPGRom]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestStructuralLSPGRom]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestFluidPGRom]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestThermalPGRom]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestStructuralPGRom]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMonotonicityPreservingRom]))

    # - testNightly
    nightlySuite = suites['nightly']
    nightlySuite.addTests(smallSuite)

    # Create a test suit that contains all the tests from every testCase
    # in the list:
    allSuite = suites['all']
    allSuite.addTests(smallSuite)

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
