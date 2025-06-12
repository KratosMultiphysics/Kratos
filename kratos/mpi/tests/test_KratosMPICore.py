# import Kratos
import KratosMultiphysics
if not KratosMultiphysics.IsDistributedRun():
    raise Exception("These tests can only be executed in MPI / distributed!")

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests or test_classes to create the suites
import test_data_communicator_factory
import test_gather_modelpart_utility
import test_mpi_communicator_set_up
import test_mpi_communicator
import test_mpi_data_communicator_python
import test_mpi_model_part
import test_mpi_serializer
import test_neighbours
import test_nodal_entity_neighbours
import test_parallel_environment
import test_mpi_processes
import test_distributed_model_part_initializer
import test_distributed_import_model_part_utility
import test_distributed_sparse_matrices

# importing OpenMP tests in MPI scope.
with KratosUnittest.WorkFolderScope("../../tests", __file__, True):
    from test_processes import TestProcesses
    import test_normal_utils
    import test_skin_detection_process
    import test_sensitivity_utilities
    from test_model_part_io import TestModelPartIOMPI
    import test_variable_redistribution
    import test_container_expression
    import test_combine_model_part_modeler
    from test_model_part_operation_utilities import TestModelPartOperationUtilities
    import test_stl_io
    import test_compute_nodal_gradient_process

def AssembleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should pupulate the suites:
    "mpi_small", "mpi_nighlty" and "mpi_all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''

    suites = KratosUnittest.KratosSuites

    # Create a test suite with the selected tests (Small tests):
    smallSuite = suites['mpi_small']
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_data_communicator_factory.TestDataCommunicatorFactory]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_gather_modelpart_utility.TestGatherModelPartUtility]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_mpi_communicator_set_up.TestMPICommunicatorSetUp]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_mpi_communicator.TestMPICommunicator]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_mpi_data_communicator_python.TestMPIDataCommunicatorPython]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_mpi_model_part.TestMPIModelPart]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_mpi_serializer.TestMPISerializer]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_neighbours.TestNeighbours]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nodal_entity_neighbours.TestNodalEntityNeighbours]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_parallel_environment.TestParallelEnvironment]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_mpi_processes.TestMPIProcesses]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_distributed_model_part_initializer.TestDistributedModelPartInitializer]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_distributed_import_model_part_utility.TestDistributedImportModelPartUtility]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_distributed_sparse_matrices.TestDistributedSparseMatrices]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_container_expression.TestHistoricalContainerExpression]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_container_expression.TestNodalContainerExpression]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_container_expression.TestConditionContainerExpression]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_container_expression.TestElementContainerExpression]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_container_expression.TestNodalPositionExpressionIO]))

    # adding non-mpi tests also as mpi tests
    smallSuite.addTest(TestProcesses("test_FindGlobalNodalNeighboursProcess"))
    smallSuite.addTest(TestProcesses("test_FindGlobalNodalNeighboursForConditionsProcess"))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_normal_utils.TestNormalUtilsCoarseSphere]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_normal_utils.TestNormalUtilsQuadSphere]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_skin_detection_process.TestSkinDetectionProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_sensitivity_utilities.TestSensitivityUtilitiesTwoDimSymmetricalSquare]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestModelPartIOMPI]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_variable_redistribution.TestVariableRedistributionUtility]))
    smallSuite.addTest(TestModelPartOperationUtilities("test_Sum"))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_combine_model_part_modeler.TestCombineModelPartModeler]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_stl_io.TestStlIO]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_compute_nodal_gradient_process.TestComputeNodalGradientProcessCoarseSphere]))

    # Create a test suite with the selected tests plus all small tests
    nightSuite = suites['mpi_nightly']
    nightSuite.addTests(smallSuite)

    # Create a test suite that contains all the tests:
    allSuite = suites['mpi_all']
    allSuite.addTests(nightSuite) # already contains the smallSuite

    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.runTests(AssembleTestSuites())
