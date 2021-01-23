from __future__ import print_function, absolute_import, division

# import Kratos
from KratosMultiphysics import *
from KratosMultiphysics.TrilinosApplication import *


# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suites

import test_trilinos_linear_solvers
import test_trilinos_matrix
import test_trilinos_redistance
import test_trilinos_levelset_convection

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

    # Create a test suite with the selected tests (Small tests):
    smallSuite = suites['mpi_small']
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_trilinos_linear_solvers.TestLinearSolvers]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_trilinos_matrix.TestTrilinosMatrix]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_trilinos_redistance.TestTrilinosRedistance]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_trilinos_levelset_convection.TestTrilinosLevelSetConvection]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_trilinos_levelset_convection.TestTrilinosLevelSetConvectionInMemory]))

    # Create a test suite with the selected tests plus all small tests
    nightSuite = suites['mpi_nightly']
    nightSuite.addTests(smallSuite)

    # Create a test suite that contains all the tests:
    allSuite = suites['mpi_all']
    allSuite.addTests(nightSuite) # already contains the smallSuite

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
