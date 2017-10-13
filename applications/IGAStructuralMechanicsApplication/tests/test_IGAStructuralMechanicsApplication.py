# import Kratos
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *
from KratosMultiphysics.IGAStructuralMechanicsApplication import *

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
## SMALL TESTS
from SmallTests import Shell_Bending_Test as TShell_Bending_Test
from SmallTests import Plate_Bending_Test as TPlate_Bending_Test
from SmallTests import Shell_Bending_Test_Lagrange as TShell_Bending_Test_Lagrange
#from SmallTests import 3PatchesTestCases as T3PatchesTestCases



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
    # smallSuite will contain the following tests:
    # - testSmallExample
    smallSuite = suites['small']
    smallSuite.addTest(TShell_Bending_Test('test_execution'))
    smallSuite.addTest(TPlate_Bending_Test('test_execution'))
    smallSuite.addTest(TShell_Bending_Test_Lagrange('test_execution'))
    #smallSuite.addTest(T3PatchesTestCases('test_execution'))

    # Create a test suit with the selected tests
    # nightSuite will contain the following tests:
    # - testSmallExample
    # - testNightlyFirstExample
    # - testNightlySecondExample
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    # Create a test suit that contains all the tests from every testCase
    # in the list:
    allSuite = suites['all']
    allSuite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            TShell_Bending_Test,
			TPlate_Bending_Test,
			TShell_Bending_Test_Lagrange
			#T3PatchesTestCases
        ])
    )

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssambleTestSuites())
