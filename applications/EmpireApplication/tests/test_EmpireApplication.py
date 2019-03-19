# import Kratos
import KratosMultiphysics
import KratosMultiphysics.EmpireApplication

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

from co_simulation_solver_test_factory import TestKratosSolver
from co_simulation_solver_test_factory import TestSDoFSolver
from co_simulation_solver_test_factory import TestMDoFSolver
from co_simulation_solver_test_factory import TestEmpireSolver
from co_simulation_test_factory import TestSmallCoSimulationCases
from co_simulation_test_factory import TestCoSimulationCases

import os
if "EMPIRE_API_LIBSO_ON_MACHINE" in os.environ:
    empire_available = True
else:
    empire_available = False # the EMPIRE environment is not set

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

    smallSuite = suites['small'] # These tests are executed by the continuous integration tool
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestKratosSolver]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestSDoFSolver]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMDoFSolver]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestEmpireSolver]))

    nightSuite = suites['nightly'] # These tests are executed in the nightly build
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestSmallCoSimulationCases]))

    nightSuite.addTests(smallSuite)

    ### Adding Validation Tests
    # For very long tests that should not be in nighly and you can use to validate
    validationSuite = suites['validation']
    validationSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCoSimulationCases]))

    # If EMPIRE is available then also add the tests involving EMPIRE
    if empire_available:
        pass








    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(nightSuite) # already contains the smallSuite
    allSuite.addTests(validationSuite)

    return suites


if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
