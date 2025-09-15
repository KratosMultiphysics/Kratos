# Importing the Kratos Library
import KratosMultiphysics as KM

if not KM.IsDistributedRun():
    raise Exception("This test script can only be executed in MPI!")

import KratosMultiphysics.KratosUnittest as KratosUnittest

from test_structural_mesh_motion_2d import TestCase as TTestCaseStructural2D
from test_structural_mesh_motion_3d import TestCase as TTestCaseStructural3D
from test_laplacian_mesh_motion_2d import TestCase as TTestCaseLaplacian2D
from test_laplacian_mesh_motion_3d import TestCase as TTestCaseLaplacian3D

# from test_ale_fluid_solver import ALEFluidSolverTest

def AssembleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should populate the suites:
    "mpi_small", "mpi_nightly" and "mpi_all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''
    suites = KratosUnittest.KratosSuites

    # Create a test suite with the selected tests (Small tests):
    smallMPISuite = suites['mpi_small']
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestCaseStructural2D]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestCaseStructural3D]))
    # TODO fix the tests/implementation!!!
    # smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestCaseLaplacian2D]))
    # smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestCaseLaplacian3D]))

    # Create a test suite with the selected tests plus all small tests
    nightlyMPISuite = suites['mpi_nightly']
    nightlyMPISuite.addTests(smallMPISuite)
    # TODO add this test!
    # nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([ALEFluidSolverTest]))

    # For very long tests that should not be in nightly and you can use to validate
    validationMPISuite = suites['mpi_validation']

    # Create a test suite that contains all the tests:
    allMPISuite = suites['mpi_all']
    allMPISuite.addTests([nightlyMPISuite]) # already contains the smallSuite
    allMPISuite.addTests([validationMPISuite])

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
