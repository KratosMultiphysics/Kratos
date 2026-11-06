import KratosMultiphysics as Kratos

if not Kratos.IsDistributedRun():
    raise Exception("This test script can only be executed in MPI!")

import KratosMultiphysics.KratosUnittest as KratosUnittest

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from test_mpi_dataset_export import TestMpiGatherHelpers
from test_mpi_dataset_export import TestMpiDatasetExport
from test_mpi_mesh_export import TestMpiGatherModelPart
from test_mpi_mesh_export import TestMpiMeshExport
from test_mpi_mesh_export import TestMpiCaeDatasetExport
from test_mpi_distributed_groups import TestMpiProcessGroups
from test_mpi_graph_partition import TestMpiHaloSubgraph
from test_mpi_graph_partition import TestMpiDataParallelTraining
from test_mpi_graph_partition import TestMpiDifferentiableHaloExchange
from test_mpi_cosim_surrogate import TestMpiDistributedSurrogateWrapper
from test_mpi_cosim_surrogate import TestMpiDistributedSurrogateCoupledLoop
from test_mpi_fsdp_checkpoint import TestMpiFsdpCheckpoint
from test_mpi_domain_parallel import TestMpiDomainParallel


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

    smallMPISuite = suites['mpi_small']
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMpiGatherHelpers]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMpiDatasetExport]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMpiGatherModelPart]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMpiMeshExport]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMpiCaeDatasetExport]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMpiProcessGroups]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMpiHaloSubgraph]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMpiDataParallelTraining]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMpiDifferentiableHaloExchange]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMpiDistributedSurrogateWrapper]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMpiDistributedSurrogateCoupledLoop]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMpiFsdpCheckpoint]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMpiDomainParallel]))

    nightlyMPISuite = suites['mpi_nightly']
    nightlyMPISuite.addTests(smallMPISuite)

    allMPISuite = suites['mpi_all']
    allMPISuite.addTests(nightlyMPISuite)  # already contains the smallMPISuite

    return suites


if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
