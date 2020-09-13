# import Kratos
from KratosMultiphysics import *
from KratosMultiphysics.IgaApplication import *
import run_cpp_unit_tests

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import Iga test factory tests
from iga_test_factory import SinglePatchTest as TSinglePatchTest
# Import patch tests
from iga_test_factory import MembraneSinglePatchFourPointSailLinearStatic as TMembraneSinglePatchFourPointSailLinearStatic
from iga_test_factory import MembraneSinglePatchFourPointSailNonLinearStatic as TMembraneSinglePatchFourPointSailNonLinearStatic
from iga_test_factory import MembraneSinglePatchFourPointSailImplicitDynamic as TMembraneSinglePatchFourPointSailImplicitDynamic
from iga_test_factory import MembraneMultiPatchFourPointSailLinearStatic as TMembraneMultiPatchFourPointSailLinearStatic
from iga_test_factory import MembraneMultiPatchFourPointSailNonLinearStatic as TMembraneMultiPatchFourPointSailNonLinearStatic
from iga_test_factory import MembraneMultiPatchFourPointSailImplicitDynamic as TMembraneMultiPatchFourPointSailImplicitDynamic
from iga_test_factory import FormfindingMembraneSinglePatchFourPointSail as TFormfindingMembraneSinglePatchFourPointSail
from iga_test_factory import FormfindingMembraneMultiPatchFourPointSail as TFormfindingMembraneMultiPatchFourPointSail

# Modelers tests
from test_modelers import TestModelers as TTestModelers

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
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TSinglePatchTest]))
    #Membrane
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TMembraneSinglePatchFourPointSailLinearStatic]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TMembraneSinglePatchFourPointSailNonLinearStatic]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TMembraneSinglePatchFourPointSailImplicitDynamic]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TMembraneMultiPatchFourPointSailLinearStatic]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TMembraneMultiPatchFourPointSailNonLinearStatic]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TMembraneMultiPatchFourPointSailImplicitDynamic]))
    #Formfinding
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TFormfindingMembraneSinglePatchFourPointSail]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TFormfindingMembraneMultiPatchFourPointSail]))

    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestModelers]))

    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    allSuite = suites['all']
    allSuite.addTests(nightSuite)

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
