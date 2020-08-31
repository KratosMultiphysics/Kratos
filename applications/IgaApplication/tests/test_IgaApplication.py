# import Kratos
from KratosMultiphysics import *
from KratosMultiphysics.IgaApplication import *
import run_cpp_unit_tests

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import Iga test factory tests
from iga_test_factory import SinglePatchTest as SinglePatchTest
# Modelers tests
from test_modelers import TestModelers as TTestModelers
from iga_test_factory import MembraneSinglePatchFourPointSailLinearStatic as MembraneSinglePatchFourPointSailLinearStatic
from iga_test_factory import MembraneSinglePatchFourPointSailNonLinearStatic as MembraneSinglePatchFourPointSailNonLinearStatic
from iga_test_factory import MembraneSinglePatchFourPointSailImplicitDynamic as MembraneSinglePatchFourPointSailImplicitDynamic

# Import the tests o test_classes to create the suits

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

    smallSuite = suites['small']
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([
        SinglePatchTest
        ]))

    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestModelers]))
    #Membrane
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([MembraneSinglePatchFourPointSailLinearStatic]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([MembraneSinglePatchFourPointSailNonLinearStatic]))

    nightSuite = suites['nightly']
    #Membrane
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([MembraneSinglePatchFourPointSailImplicitDynamic]))
    nightSuite.addTests(smallSuite)

    allSuite = suites['all']
    allSuite.addTests(nightSuite)

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
