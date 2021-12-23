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


def AssembleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should pupulate the suites:
    "small", "nighlty" and "all"

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
    smallSuite.addTest(TestRandomizedSVD('test_radomized_svd'))
    smallSuite.addTest(TestEmpiricalCubatureMethod('test_empirical_cubature_method'))


    # Create a test suit that contains all the tests from every testCase
    # in the list:
    allSuite = suites['all']
    allSuite.addTests(smallSuite)

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
