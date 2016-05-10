# import Kratos
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
from SmallTests import DynamicBossakTests as TDynamicBossakTests
from SmallTests import DynamicNewmarkTests as TDynamicNewmarkTests
from SmallTests import SprismMembranePatchTests as TSprismMembranePatchTests
from SmallTests import SprismBendingPatchTests as TSprismBendingPatchTests

def AssambleTestSuites():
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
    smallSuite = suites['small']
    smallSuite.addTest(TDynamicBossakTests('test_Bossak'))
    smallSuite.addTest(TDynamicNewmarkTests('test_Newmark'))
    smallSuite.addTest(TSprismMembranePatchTests('test_MembranePatch'))
    smallSuite.addTest(TSprismBendingPatchTests('test_BendingPatch'))

    # Create a test suit with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            TDynamicBossakTests,
            TDynamicNewmarkTests,
            TSprismMembranePatchTests,
            TSprismBendingPatchTests
        ])
    )

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssambleTestSuites())
