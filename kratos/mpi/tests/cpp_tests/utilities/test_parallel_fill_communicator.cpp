//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "mpi/testing/mpi_testing.h"
#include "containers/model.h"
#include "mpi/tests/test_utilities/mpi_cpp_test_utilities.h"
#include "mpi/utilities/parallel_fill_communicator.h"

namespace Kratos::Testing 
{

KRATOS_TEST_CASE_IN_SUITE(ParallelFillCommunicatorExecute, KratosMPICoreFastSuite)
{
    // The model part
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");
    
    // The data communicator
    const DataCommunicator& r_data_communicator = Testing::GetDefaultDataCommunicator();
    
    // MPI data
    const int rank =  r_data_communicator.Rank();
    const int world_size = r_data_communicator.Size();

    // Fill the model part
    MPICppTestUtilities::GenerateDistributedBarStructure(r_model_part, r_data_communicator);

    // Check that the communicator is correctly filled
    const auto& r_neighbours_indices = r_model_part.GetCommunicator().NeighbourIndices();
    if (world_size == 1) {
        KRATOS_EXPECT_EQ(r_neighbours_indices.size(), 0);
        KRATOS_EXPECT_EQ(r_model_part.GetCommunicator().LocalMesh().NumberOfNodes(), 11);
        KRATOS_EXPECT_EQ(r_model_part.GetCommunicator().GhostMesh().NumberOfNodes(), 0);
        KRATOS_EXPECT_EQ(r_model_part.GetCommunicator().InterfaceMesh().NumberOfNodes(), 0);
    } else if (world_size == 2) {
        KRATOS_EXPECT_EQ(r_neighbours_indices.size(), 1);
        if (rank == 0) {
            KRATOS_EXPECT_EQ(r_neighbours_indices[0], 1);
            KRATOS_EXPECT_EQ(r_model_part.GetCommunicator().LocalMesh().NumberOfNodes(), 8);
            KRATOS_EXPECT_EQ(r_model_part.GetCommunicator().GhostMesh().NumberOfNodes(), 2);
            KRATOS_EXPECT_EQ(r_model_part.GetCommunicator().InterfaceMesh().NumberOfNodes(), 2);
        } else if (rank == 1) {
            KRATOS_EXPECT_EQ(r_neighbours_indices[0], 0);
            KRATOS_EXPECT_EQ(r_model_part.GetCommunicator().LocalMesh().NumberOfNodes(), 3);
            KRATOS_EXPECT_EQ(r_model_part.GetCommunicator().GhostMesh().NumberOfNodes(), 0);
            KRATOS_EXPECT_EQ(r_model_part.GetCommunicator().InterfaceMesh().NumberOfNodes(), 2);
        }
    } else {
        if (rank == 0) {
            KRATOS_EXPECT_EQ(r_model_part.GetCommunicator().LocalMesh().NumberOfNodes(), 8);
            KRATOS_EXPECT_EQ(r_model_part.GetCommunicator().GhostMesh().NumberOfNodes(), 2);
            KRATOS_EXPECT_EQ(r_model_part.GetCommunicator().InterfaceMesh().NumberOfNodes(), 2);
        } else if (rank == 1) {
            KRATOS_EXPECT_EQ(r_model_part.GetCommunicator().LocalMesh().NumberOfNodes(), 3);
            KRATOS_EXPECT_EQ(r_model_part.GetCommunicator().GhostMesh().NumberOfNodes(), 0);
            KRATOS_EXPECT_EQ(r_model_part.GetCommunicator().InterfaceMesh().NumberOfNodes(), 2);
        } else { // The rest of the ranks
            KRATOS_EXPECT_EQ(r_model_part.GetCommunicator().LocalMesh().NumberOfNodes(), 0);
            KRATOS_EXPECT_EQ(r_model_part.GetCommunicator().GhostMesh().NumberOfNodes(), 0);
            KRATOS_EXPECT_EQ(r_model_part.GetCommunicator().InterfaceMesh().NumberOfNodes(), 0);
        }
    }
}

} // namespace Kratos::Testing