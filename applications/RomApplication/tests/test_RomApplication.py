# import Kratos
import KratosMultiphysics

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
from thermal_static_test_files.test_ROM import ROMStationaryConvDiff
from thermal_dynamic_test_files.test_ROM import ROMDynamicConvDiff
from structural_static_test_files.test_ROM import ROMStaticStruct
from structural_dynamic_test_files.test_ROM import ROMDynamicStruct
from fluid_dynamics_test_files.test_ROM import ROMFluidDynamics
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
    smallSuite.addTest(ROMStationaryConvDiff('test_ConvDiff_Stationary_ROM_2D'))
    smallSuite.addTest(ROMDynamicConvDiff('test_ConvDiff_Dynamic_ROM_2D'))
    smallSuite.addTest(ROMStaticStruct('test_Struct_Static_ROM_2D'))
    smallSuite.addTest(ROMStaticStruct('test_Struct_Static_HROM_2D'))
    smallSuite.addTest(ROMDynamicStruct('test_Struct_Dynamic_ROM_2D'))
    smallSuite.addTest(ROMFluidDynamics('test_Fluid_Dynamics_ROM_2D'))
    smallSuite.addTest(TestRandomizedSVD('test_radomized_svd'))
    smallSuite.addTest(TestEmpiricalCubatureMethod('test_empirical_cubature_method'))


    # Create a test suit that contains all the tests from every testCase
    # in the list:
    allSuite = suites['all']
    allSuite.addTests(smallSuite)

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
