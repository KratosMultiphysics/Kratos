# import Kratos
import KratosMultiphysics
import KratosMultiphysics.EigenSolversApplication as EigenSolversApplication

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
from test_eigen_direct_solver import TestEigenDirectSolver as TTestEigenDirectSolver
from test_eigensystem_solver import TestEigensystemSolver as TTestEigenSystemSolver


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
    smallSuite = suites['small']
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestEigenDirectSolver]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestEigenSystemSolver]))

    # Create a test suit with the selected tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    # Create a test suit that contains all the tests from every testCase
    # in the list:
    allSuite = suites['all']
    allSuite.addTests(nightSuite) # already contains the smallSuite

    return suites


if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
