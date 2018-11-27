# import Kratos
import KratosMultiphysics
import KratosMultiphysics.MeshMovingApplication
import run_cpp_unit_tests

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

import subprocess

from test_structural_mesh_motion_2d import TestCase as TTestCaseStructural2D
from test_structural_mesh_motion_3d import TestCase as TTestCaseStructural3D
from test_laplacian_mesh_motion_2d import TestCase as TTestCaseLaplacian2D
from test_laplacian_mesh_motion_3d import TestCase as TTestCaseLaplacian3D

## NIGTHLY TESTS

## VALIDATION TESTS

def AssembleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should populate the suites:
    "small", "nightly" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''
    suites = KratosUnittest.KratosSuites

    # Create a test suite with the selected tests (Small tests):
    smallSuite = suites['small']
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestCaseStructural2D]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestCaseStructural3D]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestCaseLaplacian2D]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestCaseLaplacian3D]))

    # Create a test suite with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    # For very long tests that should not be in nightly and you can use to validate
    validationSuite = suites['validation']

    # Create a test suite that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests([nightSuite]) # already contains the smallSuite
    allSuite.addTests([validationSuite])

    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning cpp unit tests ...")
    run_cpp_unit_tests.run()
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished running cpp unit tests!")

    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning mpi python tests ...")
    try:
        import KratosMultiphysics.mpi as KratosMPI
        import KratosMultiphysics.MetisApplication as MetisApplication
        import KratosMultiphysics.TrilinosApplication as TrilinosApplication
        p = subprocess.Popen(["mpiexec", "-np", "2", "python3", "test_MeshMovingApplication_mpi.py"], stdout=subprocess.PIPE)
        p.wait()
        KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished mpi python tests!")
    except ImportError:
        KratosMultiphysics.Logger.PrintInfo("Unittests", "mpi is not available!")

    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning python tests ...")
    KratosUnittest.runTests(AssembleTestSuites())
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished python tests!")
