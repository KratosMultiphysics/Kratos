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
import test_nodal_elemental_neighbours
import test_parallel_environment

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
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nodal_elemental_neighbours.TestNodalElementalNeighbours]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_parallel_environment.TestParallelEnvironment]))
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
