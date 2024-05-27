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
from test_affine_transform import AffineTransformTest
from test_parametric_affine_transform import ParametricAffineTransformTest
from test_impose_mesh_motion_process import TestImposeMeshMotionProcess

from test_ale_fluid_solver import ALEFluidSolverTest

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
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([AffineTransformTest]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([ParametricAffineTransformTest]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestImposeMeshMotionProcess]))

    # Create a test suite with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    # For very long tests that should not be in nightly and you can use to validate
    validationSuite = suites['validation']
    validationSuite.addTests(nightSuite)
    validationSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([ALEFluidSolverTest]))

    # Create a test suite that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests([validationSuite])

    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning cpp unit tests ...")
    run_cpp_unit_tests.run()
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished running cpp unit tests!")
    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning python tests ...")
    KratosUnittest.runTests(AssembleTestSuites())
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished python tests!")
