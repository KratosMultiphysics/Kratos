from __future__ import print_function, absolute_import, division

# import Kratos
from KratosMultiphysics import *

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests or test_classes to create the suites
import test_model_part
import test_model_part_io
import test_kratos_parameters
import test_materials_input
import test_geometries
import test_linear_solvers
import test_eigen_solvers
import test_condition_number
import test_processes
import test_importing
import test_connectivity_preserve_modeler
import test_model
import test_redistance
import test_levelset_convection
import test_variable_utils
import test_reorder
import test_exact_integration
import test_gid_io
import test_vector_interface
import test_matrix_interface
import test_restart
import test_gid_io_gauss_points
import test_skin_detection_process
import test_sparse_multiplication
import test_variable_component
import test_variable_redistribution
import test_object_printing


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
    smallSuite = suites['small']
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_model_part.TestModelPart]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_model_part_io.TestModelPartIO]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_materials_input.TestMaterialsInput]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_geometries.TestGeometry]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_kratos_parameters.TestParameters]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_linear_solvers.TestLinearSolvers]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_eigen_solvers.TestEigenSolvers]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_condition_number.TestConditionNumber]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_processes.TestProcesses]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_importing.TestImporting]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_connectivity_preserve_modeler.TestConnectivityPreserveModeler]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_model.TestModel]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_redistance.TestRedistance]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_levelset_convection.TestLevelSetConvection]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_variable_utils.TestVariableUtils]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_reorder.TestReorder]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_exact_integration.TestExactIntegration]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_gid_io.TestGidIO]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_vector_interface.TestVectorInterface]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_matrix_interface.TestMatrixInterface]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_restart.TestRestart]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_gid_io_gauss_points.TestGiDIOGaussPoints]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_skin_detection_process.TestSkinDetectionProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_sparse_multiplication.TestSparseMatrixSum]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_sparse_multiplication.TestSparseMatrixTranspose]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_sparse_multiplication.TestSparseMatrixMultiplication]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_variable_component.TestVariableComponent]))
    smallSuite.addTest(test_variable_redistribution.VariableRedistributionTest('testLinearFunction'))
    smallSuite.addTest(test_variable_redistribution.VariableRedistributionTest('testSharpCorners'))
    smallSuite.addTest(test_variable_redistribution.VariableRedistributionTest('testVector'))
    smallSuite.addTest(test_variable_redistribution.VariableRedistributionTest('testQuadratic'))
    smallSuite.addTest(test_variable_redistribution.VariableRedistributionTest('testNodalArea'))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_object_printing.TestObjectPrinting]))

    # Create a test suite with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    # Create a test suite that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(nightSuite) # already contains the smallSuite

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
