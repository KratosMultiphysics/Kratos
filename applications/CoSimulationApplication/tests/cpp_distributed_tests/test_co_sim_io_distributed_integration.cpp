// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:		 BSD License
//                   license: CoSimulationApplication/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/parallel_environment.h"
#include "containers/model.h"
#include "tests/cpp_distributed_tests/co_simulation_distributed_suite.h"
#include "custom_utilities/co_sim_io_conversion_utilities.h"
#include "tests/test_utilities/co_simulation_testing_utilities.h"

// Utilities

namespace Kratos {
namespace Testing {

using DataLocation = Globals::DataLocation;

namespace DistributedTestHelpers {

std::size_t GetPartnerRank()
{
    const auto& r_world_data_comm = ParallelEnvironment::GetDataCommunicator("World");
    const std::size_t my_rank = r_world_data_comm.Rank();
    const std::size_t world_size = r_world_data_comm.Size();

    return (my_rank+1) % world_size;
}

std::size_t GetId(const std::size_t NumLocalNodesPerRank, const std::size_t LocalNodeIndex)
{
    const std::size_t my_rank = ParallelEnvironment::GetDataCommunicator("World").Rank();
    return my_rank*NumLocalNodesPerRank + LocalNodeIndex+1;
}

std::size_t GetGhostId(const std::size_t NumLocalNodesPerRank, const std::size_t LocalNodeIndex)
{
    return GetPartnerRank()*NumLocalNodesPerRank + LocalNodeIndex+1;
}

void CreateDistributedNodes(
    CoSimIO::ModelPart& rCoSimIOModelPart,
    const std::size_t NumLocalNodesPerRank,
    const std::size_t NumGhostNodesPerRank)
{
    KRATOS_EXPECT_GT(NumLocalNodesPerRank, NumGhostNodesPerRank);

    const auto& r_world_data_comm = ParallelEnvironment::GetDataCommunicator("World");
    const std::size_t my_rank = r_world_data_comm.Rank();
    const std::size_t world_size = r_world_data_comm.Size();

    auto create_ghost_nodes = [&](){
        for (std::size_t i=0; i<NumGhostNodesPerRank; ++i) {
            rCoSimIOModelPart.CreateNewGhostNode(GetGhostId(NumLocalNodesPerRank, i), 0,0,0, GetPartnerRank());
        }
    };

    if (GetPartnerRank() < my_rank) {
        // Kratos requires a consecutively ordered input as it sorts the nodes internally and then the matching would be lost
        // in case the ghost nodes have smaller ids than the local nodes they have to be added first as currently
        // the CoSimIO::ModelPart does not allow sorting
        create_ghost_nodes();
    }

    // create local nodes
    for (std::size_t i=0; i<NumLocalNodesPerRank; ++i) {
        rCoSimIOModelPart.CreateNewNode(GetId(NumLocalNodesPerRank, i), 0,0,0);
    }

    if (GetPartnerRank() > my_rank) {
        // in this case the ghost nodes have larger Ids than the local nodes
        create_ghost_nodes();
    }

    KRATOS_EXPECT_EQ(rCoSimIOModelPart.NumberOfNodes(), NumLocalNodesPerRank+NumGhostNodesPerRank);
    KRATOS_EXPECT_EQ(rCoSimIOModelPart.NumberOfLocalNodes(), NumLocalNodesPerRank);
    KRATOS_EXPECT_EQ(rCoSimIOModelPart.NumberOfGhostNodes(), NumGhostNodesPerRank);
    KRATOS_EXPECT_EQ(r_world_data_comm.SumAll(static_cast<int>(rCoSimIOModelPart.NumberOfLocalNodes())), static_cast<int>(NumLocalNodesPerRank*world_size));
    KRATOS_EXPECT_EQ(rCoSimIOModelPart.NumberOfElements(), 0);
}

void CreateDistributedNodesUnordered(
    CoSimIO::ModelPart& rCoSimIOModelPart,
    const std::size_t NumLocalNodesPerRank,
    const std::size_t NumGhostNodesPerRank)
{
    KRATOS_EXPECT_GT(NumLocalNodesPerRank, NumGhostNodesPerRank);

    const auto& r_world_data_comm = ParallelEnvironment::GetDataCommunicator("World");
    const std::size_t world_size = r_world_data_comm.Size();

    for (std::size_t i=0; i<NumLocalNodesPerRank; ++i) {
        rCoSimIOModelPart.CreateNewNode(GetId(NumLocalNodesPerRank, i), 0,0,0);
    }

    for (int i=static_cast<int>(NumGhostNodesPerRank); --i>=0;) { // must use signed type here!
        rCoSimIOModelPart.CreateNewGhostNode(GetGhostId(NumLocalNodesPerRank, i), 0,0,0, GetPartnerRank());
    }

    KRATOS_EXPECT_EQ(rCoSimIOModelPart.NumberOfNodes(), NumLocalNodesPerRank+NumGhostNodesPerRank);
    KRATOS_EXPECT_EQ(rCoSimIOModelPart.NumberOfLocalNodes(), NumLocalNodesPerRank);
    KRATOS_EXPECT_EQ(rCoSimIOModelPart.NumberOfGhostNodes(), NumGhostNodesPerRank);
    KRATOS_EXPECT_EQ(r_world_data_comm.SumAll(static_cast<int>(rCoSimIOModelPart.NumberOfLocalNodes())), static_cast<int>(NumLocalNodesPerRank*world_size));
    KRATOS_EXPECT_EQ(rCoSimIOModelPart.NumberOfElements(), 0);
}

void CreateDistributedNodes(
    Kratos::ModelPart& rKratosModelPart,
    const std::size_t NumLocalNodesPerRank,
    const std::size_t NumGhostNodesPerRank)
{
    KRATOS_EXPECT_GT(NumLocalNodesPerRank, NumGhostNodesPerRank);

    const auto& r_world_data_comm = ParallelEnvironment::GetDataCommunicator("World");
    const std::size_t my_rank = r_world_data_comm.Rank();
    const std::size_t world_size = r_world_data_comm.Size();

    // create local nodes
    for (std::size_t i=0; i<NumLocalNodesPerRank; ++i) {
        auto node = rKratosModelPart.CreateNewNode(GetId(NumLocalNodesPerRank, i), i*10,i/3.58,i-99.14);
        node->FastGetSolutionStepValue(PARTITION_INDEX) = my_rank;
    }

    // create ghost nodes
    for (std::size_t i=0; i<NumGhostNodesPerRank; ++i) {
        auto node = rKratosModelPart.CreateNewNode(GetGhostId(NumLocalNodesPerRank, i), i-10,i*3,i+99.14);
        node->FastGetSolutionStepValue(PARTITION_INDEX) = GetPartnerRank();
    }

    // this calls the ParallelFillCommunicator
    ParallelEnvironment::CreateFillCommunicatorFromGlobalParallelism(rKratosModelPart, r_world_data_comm)->Execute();

    KRATOS_EXPECT_EQ(rKratosModelPart.NumberOfNodes(), NumLocalNodesPerRank+NumGhostNodesPerRank);
    KRATOS_EXPECT_EQ(rKratosModelPart.GetCommunicator().LocalMesh().NumberOfNodes(), NumLocalNodesPerRank);
    KRATOS_EXPECT_EQ(rKratosModelPart.GetCommunicator().GhostMesh().NumberOfNodes(), NumGhostNodesPerRank);
    KRATOS_EXPECT_EQ(rKratosModelPart.GetCommunicator().GlobalNumberOfNodes(), (NumLocalNodesPerRank)*world_size);
    KRATOS_EXPECT_EQ(rKratosModelPart.NumberOfProperties(), 0);
}

} // DistributedTestHelpers

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(DistributedCoSimIOModelPartToKratosModelPart_NodesOnly, KratosCoSimulationMPIFastSuite)
{
    const auto& r_world_data_comm = ParallelEnvironment::GetDataCommunicator("World");

    Model model;
    auto& kratos_model_part = model.CreateModelPart("kratos_mp");

    CoSimIO::ModelPart co_sim_io_model_part("co_sim_io_mp");

    constexpr std::size_t num_local_nodes_per_rank = 5;
    constexpr std::size_t num_ghost_nodes_per_rank = 3;

    DistributedTestHelpers::CreateDistributedNodes(co_sim_io_model_part, num_local_nodes_per_rank, num_ghost_nodes_per_rank);

    CoSimIOConversionUtilities::CoSimIOModelPartToKratosModelPart(co_sim_io_model_part, kratos_model_part, r_world_data_comm);

    CheckDistributedModelPartsAreEqual(kratos_model_part, co_sim_io_model_part);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(DistributedCoSimIOModelPartToKratosModelPart_NodesOnly_Unordered, KratosCoSimulationMPIFastSuite)
{
    const auto& r_world_data_comm = ParallelEnvironment::GetDataCommunicator("World");

    Model model;
    auto& kratos_model_part = model.CreateModelPart("kratos_mp");

    CoSimIO::ModelPart co_sim_io_model_part("co_sim_io_mp");

    constexpr std::size_t num_local_nodes_per_rank = 5;
    constexpr std::size_t num_ghost_nodes_per_rank = 3;

    DistributedTestHelpers::CreateDistributedNodesUnordered(co_sim_io_model_part, num_local_nodes_per_rank, num_ghost_nodes_per_rank);

    CoSimIOConversionUtilities::CoSimIOModelPartToKratosModelPart(co_sim_io_model_part, kratos_model_part, r_world_data_comm);

    CheckDistributedModelPartsAreEqualButEntitiesAreOrderedDifferently(kratos_model_part, co_sim_io_model_part);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(DistributedKratosModelPartToCoSimIOModelPart_NodesOnly, KratosCoSimulationMPIFastSuite)
{
    Model model;
    auto& kratos_model_part = model.CreateModelPart("kratos_mp");
    kratos_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    CoSimIO::ModelPart co_sim_io_model_part("co_sim_io_mp");

    constexpr std::size_t num_local_nodes_per_rank = 5;
    constexpr std::size_t num_ghost_nodes_per_rank = 3;

    DistributedTestHelpers::CreateDistributedNodes(kratos_model_part, num_local_nodes_per_rank, num_ghost_nodes_per_rank);

    CoSimIOConversionUtilities::KratosModelPartToCoSimIOModelPart(kratos_model_part, co_sim_io_model_part);

    CheckDistributedModelPartsAreEqual(kratos_model_part, co_sim_io_model_part);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(DistributedCoSimIOModelPartToKratosModelPart, KratosCoSimulationMPIFastSuite)
{
    const auto& r_world_data_comm = ParallelEnvironment::GetDataCommunicator("World");

    Model model;
    auto& kratos_model_part = model.CreateModelPart("kratos_mp");

    CoSimIO::ModelPart co_sim_io_model_part("co_sim_io_mp");

    constexpr std::size_t num_local_nodes_per_rank = 5;
    constexpr std::size_t num_ghost_nodes_per_rank = 3;

    DistributedTestHelpers::CreateDistributedNodes(co_sim_io_model_part, num_local_nodes_per_rank, num_ghost_nodes_per_rank);

    // elements that use two local nodes
    for (std::size_t i=0; i<num_ghost_nodes_per_rank; ++i) {
        CoSimIO::ConnectivitiesType conn {
            static_cast<int>(DistributedTestHelpers::GetId(num_local_nodes_per_rank, i)),
            static_cast<int>(DistributedTestHelpers::GetId(num_local_nodes_per_rank, i+1))};
        co_sim_io_model_part.CreateNewElement(DistributedTestHelpers::GetId(num_local_nodes_per_rank, i), CoSimIO::ElementType::Line2D2, conn);
    }

    // elements that use one local and one ghost node
    for (std::size_t i=0; i<num_ghost_nodes_per_rank; ++i) {
        CoSimIO::ConnectivitiesType conn {
            static_cast<int>(DistributedTestHelpers::GetId(num_local_nodes_per_rank, i)),
            static_cast<int>(DistributedTestHelpers::GetGhostId(num_local_nodes_per_rank, i))};
        co_sim_io_model_part.CreateNewElement(DistributedTestHelpers::GetId(num_local_nodes_per_rank, i)+num_ghost_nodes_per_rank, CoSimIO::ElementType::Line2D2, conn);
    }

    CoSimIOConversionUtilities::CoSimIOModelPartToKratosModelPart(co_sim_io_model_part, kratos_model_part, r_world_data_comm);

    CheckDistributedModelPartsAreEqual(kratos_model_part, co_sim_io_model_part);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(DistributedCoSimIOModelPartToKratosModelPart_Unordered, KratosCoSimulationMPIFastSuite)
{
    const auto& r_world_data_comm = ParallelEnvironment::GetDataCommunicator("World");

    Model model;
    auto& kratos_model_part = model.CreateModelPart("kratos_mp");

    CoSimIO::ModelPart co_sim_io_model_part("co_sim_io_mp");

    constexpr std::size_t num_local_nodes_per_rank = 5;
    constexpr std::size_t num_ghost_nodes_per_rank = 3;

    DistributedTestHelpers::CreateDistributedNodesUnordered(co_sim_io_model_part, num_local_nodes_per_rank, num_ghost_nodes_per_rank);

    // elements that use two local nodes
    for (std::size_t i=0; i<num_ghost_nodes_per_rank; ++i) {
        CoSimIO::ConnectivitiesType conn {
            static_cast<int>(DistributedTestHelpers::GetId(num_local_nodes_per_rank, i)),
            static_cast<int>(DistributedTestHelpers::GetId(num_local_nodes_per_rank, i+1))};
        co_sim_io_model_part.CreateNewElement(DistributedTestHelpers::GetId(num_local_nodes_per_rank, i), CoSimIO::ElementType::Line2D2, conn);
    }

    // elements that use one local and one ghost node
    for (std::size_t i=0; i<num_ghost_nodes_per_rank; ++i) {
        CoSimIO::ConnectivitiesType conn {
            static_cast<int>(DistributedTestHelpers::GetId(num_local_nodes_per_rank, i)),
            static_cast<int>(DistributedTestHelpers::GetGhostId(num_local_nodes_per_rank, i))};
        co_sim_io_model_part.CreateNewElement(DistributedTestHelpers::GetId(num_local_nodes_per_rank, i)+num_ghost_nodes_per_rank, CoSimIO::ElementType::Line2D2, conn);
    }

    CoSimIOConversionUtilities::CoSimIOModelPartToKratosModelPart(co_sim_io_model_part, kratos_model_part, r_world_data_comm);

    CheckDistributedModelPartsAreEqualButEntitiesAreOrderedDifferently(kratos_model_part, co_sim_io_model_part);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(KratosDistributedModelPartToCoSimIOModelPart, KratosCoSimulationMPIFastSuite)
{
    const auto& r_world_data_comm = ParallelEnvironment::GetDataCommunicator("World");

    Model model;
    auto& kratos_model_part = model.CreateModelPart("kratos_mp");
    kratos_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    CoSimIO::ModelPart co_sim_io_model_part("co_sim_io_mp");

    constexpr std::size_t num_local_nodes_per_rank = 5;
    constexpr std::size_t num_ghost_nodes_per_rank = 3;

    DistributedTestHelpers::CreateDistributedNodes(kratos_model_part, num_local_nodes_per_rank, num_ghost_nodes_per_rank);

    Properties::Pointer p_props(kratos_model_part.CreateNewProperties(0));

    // elements that use two local nodes
    for (std::size_t i=0; i<num_ghost_nodes_per_rank; ++i) {
        std::vector<ModelPart::IndexType> conn {
            static_cast<ModelPart::IndexType>(DistributedTestHelpers::GetId(num_local_nodes_per_rank, i)),
            static_cast<ModelPart::IndexType>(DistributedTestHelpers::GetId(num_local_nodes_per_rank, i+1))};
        const std::size_t elem_id = DistributedTestHelpers::GetId(num_local_nodes_per_rank, i);
        kratos_model_part.CreateNewElement("Element2D2N", elem_id, conn, p_props);
    }

    // elements that use one local and one ghost node
    for (std::size_t i=0; i<num_ghost_nodes_per_rank; ++i) {
        std::vector<ModelPart::IndexType> conn {
            static_cast<ModelPart::IndexType>(DistributedTestHelpers::GetId(num_local_nodes_per_rank, i)),
            static_cast<ModelPart::IndexType>(DistributedTestHelpers::GetGhostId(num_local_nodes_per_rank, i))};
        const std::size_t elem_id = DistributedTestHelpers::GetId(num_local_nodes_per_rank, i)+num_ghost_nodes_per_rank;
        kratos_model_part.CreateNewElement("Element2D2N", elem_id, conn, p_props);
    }

    // this calls the ParallelFillCommunicator
    ParallelEnvironment::CreateFillCommunicatorFromGlobalParallelism(kratos_model_part, r_world_data_comm)->Execute();

    CoSimIOConversionUtilities::KratosModelPartToCoSimIOModelPart(kratos_model_part, co_sim_io_model_part);

    CheckDistributedModelPartsAreEqual(kratos_model_part, co_sim_io_model_part);
}

} // namespace Testing
} // namespace Kratos
