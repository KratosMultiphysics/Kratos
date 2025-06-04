//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi

// System includes

// External includes

// Project includes
#include "mpi/testing/mpi_testing.h"
#include "mpi/tests/test_utilities/distributed_sparse_containers_test_utilities.h"
#include "containers/sparse_contiguous_row_graph.h"
#include "containers/distributed_sparse_graph.h"
#include "utilities/builtin_timer.h"

#include "mpi/utilities/amgcl_distributed_csr_conversion_utilities.h"
#include "mpi/utilities/amgcl_distributed_csr_spmm_utilities.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(DistributedSparseGraphConstruction, KratosMPICoreFastSuite)
{
    using IndexType = DistributedSparseContainersTestUtilities::IndexType;
    const auto& r_comm = ParallelEnvironment::GetDefaultDataCommunicator();
    int world_size = r_comm.Size();
    int my_rank = r_comm.Rank();

    auto dofs_bounds = DistributedSparseContainersTestUtilities::ComputeBounds<IndexType>(40, world_size, my_rank);
    auto reference_A_map = DistributedSparseContainersTestUtilities::GetReferenceMatrixAsMap(dofs_bounds);

    auto el_bounds = DistributedSparseContainersTestUtilities::ComputeBounds<IndexType>(31, world_size, my_rank);
    const auto connectivities = DistributedSparseContainersTestUtilities::ElementConnectivities(el_bounds);

    DistributedSparseGraph<IndexType> Agraph(dofs_bounds[1]-dofs_bounds[0], r_comm);

    IndexPartition<IndexType>(connectivities.size()).for_each([&](DistributedSparseContainersTestUtilities::IndexType i) {
        Agraph.AddEntries(connectivities[i]);
    });
    Agraph.Finalize();

    DistributedSparseContainersTestUtilities::CheckGraph(Agraph, reference_A_map);
}

KRATOS_TEST_CASE_IN_SUITE(DistributedSparseGraphConstructionBenchmark, KratosMPICoreFastSuite)
{
    using IndexType = DistributedSparseContainersTestUtilities::IndexType;
    const auto& r_comm = ParallelEnvironment::GetDefaultDataCommunicator();
    int world_size = r_comm.Size();
    int my_rank = r_comm.Rank();

    const IndexType block_size = 4;
    const IndexType nodes_in_elem = 4;
    const IndexType nel = 1e2; //set to 1e6 or 1e7 for a more realistic benchmark
    const IndexType ndof = nel/6;
    const IndexType standard_dev = 100; //reducing this implies using more optimized numbering

    auto el_bounds = DistributedSparseContainersTestUtilities::ComputeBounds<DistributedSparseContainersTestUtilities::IndexType>(nel, world_size, my_rank);
    auto dofs_bounds = DistributedSparseContainersTestUtilities::ComputeBounds<DistributedSparseContainersTestUtilities::IndexType>(ndof*block_size, world_size, my_rank);

    auto connectivities = SparseContainersTestUtilities::RandomElementConnectivities(block_size,
                          nodes_in_elem,
                          el_bounds[0],
                          el_bounds[1],
                          ndof,
                          standard_dev);

    r_comm.Barrier(); //to ensure fair timings
    // auto timer =  BuiltinTimer();
    DistributedSparseGraph<IndexType> Agraph(dofs_bounds[1]-dofs_bounds[0], r_comm);

    IndexPartition<IndexType>(connectivities.size()).for_each([&](IndexType i) {
        Agraph.AddEntries(connectivities[i]);
    });
    Agraph.Finalize();
    r_comm.Barrier(); //to ensure fair timings

    // std::cout << "DistributedSparseGraphConstructionBenchmark - time = " << timer.ElapsedSeconds() << std::endl;
}

} // namespace Kratos::Testing
