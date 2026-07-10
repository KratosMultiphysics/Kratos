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
#include "containers/distributed_system_vector.h"
#include "containers/distributed_vector_importer.h"
#include "mpi/includes/mpi_data_communicator.h"
#include "mpi/testing/mpi_testing.h"
#include "mpi/tests/test_utilities/distributed_sparse_containers_test_utilities.h"
#include "mpi/utilities/amgcl_distributed_csr_conversion_utilities.h"
#include "mpi/utilities/amgcl_distributed_csr_spmm_utilities.h"
#include "utilities/openmp_utils.h"
#include "utilities/builtin_timer.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(DistributedSystemVectorAssembly, KratosMPICoreFastSuite)
{
    const auto& r_comm = ParallelEnvironment::GetDefaultDataCommunicator();
    int world_size = r_comm.Size();
    int my_rank = r_comm.Rank();

    auto dofs_bounds = DistributedSparseContainersTestUtilities::ComputeBounds<DistributedSparseContainersTestUtilities::IndexType>(40, world_size, my_rank);
    auto reference_b_map = DistributedSparseContainersTestUtilities::GetReferenceVectorAsMap(dofs_bounds);

    auto el_bounds = DistributedSparseContainersTestUtilities::ComputeBounds<DistributedSparseContainersTestUtilities::IndexType>(31, world_size, my_rank);
    const auto connectivities = DistributedSparseContainersTestUtilities::ElementConnectivities(el_bounds);

    DistributedSparseGraph<DistributedSparseContainersTestUtilities::IndexType> dist_sparse_graph(dofs_bounds[1] - dofs_bounds[0], r_comm);
    IndexPartition<DistributedSparseContainersTestUtilities::IndexType>(connectivities.size()).for_each([&](DistributedSparseContainersTestUtilities::IndexType i) {
        dist_sparse_graph.AddEntries(connectivities[i]);
    });
    dist_sparse_graph.Finalize();

    //FEM assembly
    DistributedSystemVector<> b(dist_sparse_graph);
    b.SetValue(0.0);
    b.BeginAssemble();
    for(const auto& c : connectivities) {
        Vector vdata(c.size(),1.0);
        b.Assemble(vdata,c);
    }
    b.FinalizeAssemble();

    DistributedSparseContainersTestUtilities::IndexType local_size = b.LocalSize();
    for(unsigned int i = 0; i < local_size; ++i) {
        auto global_i = b.GetNumbering().GlobalId(i);
        auto it = reference_b_map.find(global_i);
        const auto& ref_value = it->second;
        KRATOS_EXPECT_NEAR(b(i),  ref_value, 1e-14 );
    }

    //test importing
    std::vector<DistributedSparseContainersTestUtilities::IndexType> to_import{39, 0,37,2};
    DistributedVectorImporter<double, DistributedSparseContainersTestUtilities::IndexType> importer(b.GetComm(), to_import, b.GetNumbering()); //this operation is expensive since it requires mounting the communication plan
    auto x = importer.ImportData(b);
    KRATOS_EXPECT_NEAR(x[0],  3, 1e-14 );
    KRATOS_EXPECT_NEAR(x[1],  1, 1e-14 );
    KRATOS_EXPECT_NEAR(x[2],  4, 1e-14 );
    KRATOS_EXPECT_NEAR(x[3],  2, 1e-14 );
}

KRATOS_TEST_CASE_IN_SUITE(DistributedSystemVectorOperations, KratosMPICoreFastSuite)
{
    const auto& r_comm = ParallelEnvironment::GetDefaultDataCommunicator();

    DistributedSparseContainersTestUtilities::IndexType local_size = 4;
    DistributedNumbering<DistributedSparseContainersTestUtilities::IndexType> numbering(r_comm, local_size);
    DistributedSystemVector<double,DistributedSparseContainersTestUtilities::IndexType> a(numbering);
    KRATOS_EXPECT_EQ(a.LocalSize(), local_size);
    a.SetValue(5.0);
    for(unsigned int i=0; i<a.LocalSize(); ++i) {
        KRATOS_EXPECT_NEAR(a[i], 5.0,1e-14);
    }

    DistributedSystemVector<double,DistributedSparseContainersTestUtilities::IndexType> b(numbering);
    b.SetValue(3.0);
    KRATOS_EXPECT_EQ(b.LocalSize(), local_size);
    for(unsigned int i=0; i<b.LocalSize(); ++i) {
        KRATOS_EXPECT_NEAR(b[i], 3.0,1e-14);
    }

    DistributedSystemVector<double,DistributedSparseContainersTestUtilities::IndexType> c(a);
    KRATOS_EXPECT_EQ(c.LocalSize(), local_size);
    KRATOS_EXPECT_EQ(c.Size(), local_size*r_comm.Size());
    for(unsigned int i=0; i<c.LocalSize(); ++i) {
        KRATOS_EXPECT_NEAR(c[i], 5.0,1e-14);
    }

    c += b;
    for(unsigned int i=0; i<c.LocalSize(); ++i) {
        KRATOS_EXPECT_NEAR(c[i], 8.0,1e-14);
    }

    c -= b;
    for(unsigned int i=0; i<c.LocalSize(); ++i)
        KRATOS_EXPECT_NEAR(c[i], 5.0,1e-14);

    c.Add(3.0,a);
    for(unsigned int i=0; i<c.LocalSize(); ++i) {
        KRATOS_EXPECT_NEAR(c[i], 20.0,1e-14);
    }

    c*=2.0;
    for(unsigned int i=0; i<c.LocalSize(); ++i) {
        KRATOS_EXPECT_NEAR(c[i], 40.0,1e-14);
    }

    c/=4.0;
    for(unsigned int i=0; i<c.LocalSize(); ++i) {
        KRATOS_EXPECT_NEAR(c[i], 10.0,1e-14);
    }
}

} // namespace Kratos::Testing
