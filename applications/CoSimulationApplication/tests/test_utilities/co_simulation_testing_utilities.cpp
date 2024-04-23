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
#include "co_simulation_application_variables.h"
#include "tests/test_utilities/co_simulation_testing_utilities.h"


namespace Kratos {
namespace Testing {

namespace { // helpers namespace

void CheckEntitiesAreEqual(
    const Kratos::Node& rKratosNode,
    const CoSimIO::Node& rCoSimIONode)
{
    KRATOS_TRY

    KRATOS_EXPECT_EQ(rKratosNode.Id(), static_cast<std::size_t>(rCoSimIONode.Id()));

    KRATOS_EXPECT_DOUBLE_EQ(rKratosNode.X(),  rCoSimIONode.X());
    KRATOS_EXPECT_DOUBLE_EQ(rKratosNode.X0(), rCoSimIONode.X());

    KRATOS_EXPECT_DOUBLE_EQ(rKratosNode.Y(),  rCoSimIONode.Y());
    KRATOS_EXPECT_DOUBLE_EQ(rKratosNode.Y0(), rCoSimIONode.Y());

    KRATOS_EXPECT_DOUBLE_EQ(rKratosNode.Z(),  rCoSimIONode.Z());
    KRATOS_EXPECT_DOUBLE_EQ(rKratosNode.Z0(), rCoSimIONode.Z());

    KRATOS_CATCH("")
}

void CheckEntitiesAreEqual(
    const Kratos::Element& rKratosElement,
    const CoSimIO::Element& rCoSimIOElement)
{
    KRATOS_TRY

    // basic checks
    KRATOS_EXPECT_EQ(rKratosElement.Id(), static_cast<std::size_t>(rCoSimIOElement.Id()));
    // KRATOS_EXPECT_EQ(rKratosElement.Type(), rCoSimIOElement.Type());
    KRATOS_EXPECT_EQ(rKratosElement.GetGeometry().PointsNumber(), rCoSimIOElement.NumberOfNodes());

    // check nodes
    for (std::size_t i=0; i<rCoSimIOElement.NumberOfNodes(); ++i) {
        CheckNodesAreEqual(*(rKratosElement.GetGeometry().begin()+i), **(rCoSimIOElement.NodesBegin()+i));
    }

    KRATOS_CATCH("")
}

template<class TKratosContainerType, class TCoSimIOContainerType>
void CheckUnorderedEntitiesAreEqual(
    const TKratosContainerType& rKratosEntities,
    const TCoSimIOContainerType& rCoSimIOEntities)
{
    KRATOS_TRY

    // basic checks
    KRATOS_EXPECT_EQ(rKratosEntities.size(), rCoSimIOEntities.size());

    // check entities
    for (const auto& r_cosimio_entity : rCoSimIOEntities) {
        const auto it_kratos_entity = rKratosEntities.find(r_cosimio_entity.Id());
        KRATOS_ERROR_IF(it_kratos_entity == rKratosEntities.end()) << "Entity with Id " << r_cosimio_entity.Id() << " is not in the Kratos ModelPart!" << std::endl;
        const auto& r_kratos_entity = *it_kratos_entity;
        CheckEntitiesAreEqual(r_kratos_entity, r_cosimio_entity);
    }

    KRATOS_CATCH("")
}

void CheckNumberOfGhostNodes(
    const Kratos::ModelPart& rKratosModelPart,
    const CoSimIO::ModelPart& rCoSimIOModelPart)
{
    KRATOS_TRY

    std::vector<std::size_t> kratos_ghost_nodes_per_rank(rKratosModelPart.GetCommunicator().TotalProcesses(), 0);
    for (const auto& r_kratos_node : rKratosModelPart.GetCommunicator().GhostMesh().Nodes()) {
        kratos_ghost_nodes_per_rank[r_kratos_node.FastGetSolutionStepValue(PARTITION_INDEX)]++;
    }

    for (int part_index=0; part_index<rKratosModelPart.GetCommunicator().TotalProcesses(); ++part_index) {
        const std::size_t num_kratos_nodes_this_rank = kratos_ghost_nodes_per_rank[part_index];
        if (num_kratos_nodes_this_rank > 0) {
            KRATOS_EXPECT_EQ(rCoSimIOModelPart.GetPartitionModelParts().count(part_index), 1);
            CoSimIO::ModelPart& r_cosimio_part_mp = *rCoSimIOModelPart.GetPartitionModelParts().at(part_index);
            KRATOS_EXPECT_EQ(num_kratos_nodes_this_rank, r_cosimio_part_mp.NumberOfNodes());
        } else {
            KRATOS_EXPECT_EQ(rCoSimIOModelPart.GetPartitionModelParts().count(part_index), 0);
        }
    }

    KRATOS_CATCH("")
}

} // helpers namespace

void CheckNodesAreEqual(
    const Kratos::Node& rKratosNode,
    const CoSimIO::Node& rCoSimIONode)
{
    KRATOS_TRY

    CheckEntitiesAreEqual(rKratosNode, rCoSimIONode);

    KRATOS_CATCH("")
}

void CheckElementsAreEqual(
    const Kratos::Element& rKratosElement,
    const CoSimIO::Element& rCoSimIOElement)
{
    KRATOS_TRY

    CheckEntitiesAreEqual(rKratosElement, rCoSimIOElement);

    KRATOS_CATCH("")
}

void CheckModelPartsAreEqual(
    const Kratos::ModelPart& rKratosModelPart,
    const CoSimIO::ModelPart& rCoSimIOModelPart)
{
    KRATOS_TRY

    // basic checks
    KRATOS_EXPECT_EQ(rKratosModelPart.NumberOfNodes(),      rCoSimIOModelPart.NumberOfNodes());
    KRATOS_EXPECT_EQ(rKratosModelPart.NumberOfElements(),   rCoSimIOModelPart.NumberOfElements());

    // check nodes
    for (std::size_t i=0; i<rCoSimIOModelPart.NumberOfNodes(); ++i) {
        CheckNodesAreEqual(*(rKratosModelPart.NodesBegin()+i), **(rCoSimIOModelPart.NodesBegin()+i));
    }

    // check elements
    for (std::size_t i=0; i<rCoSimIOModelPart.NumberOfElements(); ++i) {
        CheckElementsAreEqual(*(rKratosModelPart.ElementsBegin()+i), **(rCoSimIOModelPart.ElementsBegin()+i));
    }

    KRATOS_CATCH("")
}

void CheckModelPartsAreEqualButEntitiesAreOrderedDifferently(
    const Kratos::ModelPart& rKratosModelPart,
    const CoSimIO::ModelPart& rCoSimIOModelPart)
{
    KRATOS_TRY

    KRATOS_EXPECT_TRUE(rKratosModelPart.Has(NODES_ID_INDEX_MAP));
    KRATOS_EXPECT_TRUE(rKratosModelPart.Has(ELEMENTS_ID_INDEX_MAP));

    // basic checks
    KRATOS_EXPECT_EQ(rKratosModelPart.NumberOfNodes(),      rCoSimIOModelPart.NumberOfNodes());
    KRATOS_EXPECT_EQ(rKratosModelPart.NumberOfElements(),   rCoSimIOModelPart.NumberOfElements());

    const auto& r_nodes_id_map = rKratosModelPart[NODES_ID_INDEX_MAP];
    const auto& r_elements_id_map = rKratosModelPart[ELEMENTS_ID_INDEX_MAP];

    KRATOS_EXPECT_EQ(rKratosModelPart.NumberOfNodes(),    r_nodes_id_map.size());
    KRATOS_EXPECT_EQ(rKratosModelPart.NumberOfElements(), r_elements_id_map.size());

    // check nodes
    CheckUnorderedEntitiesAreEqual(rKratosModelPart.Nodes(), rCoSimIOModelPart.Nodes());

    // check elements
    CheckUnorderedEntitiesAreEqual(rKratosModelPart.Elements(), rCoSimIOModelPart.Elements());

    // check that orderings are created correctly
    for (std::size_t i=0; i<rCoSimIOModelPart.NumberOfNodes(); ++i) {
        KRATOS_EXPECT_EQ(r_nodes_id_map[i], static_cast<std::size_t>((**(rCoSimIOModelPart.NodesBegin()+i)).Id()));
    }
    for (std::size_t i=0; i<rCoSimIOModelPart.NumberOfElements(); ++i) {
        KRATOS_EXPECT_EQ(r_elements_id_map[i], static_cast<std::size_t>((**(rCoSimIOModelPart.ElementsBegin()+i)).Id()));
    }

    KRATOS_CATCH("")
}

void CheckDistributedModelPartsAreEqual(
    const Kratos::ModelPart& rKratosModelPart,
    const CoSimIO::ModelPart& rCoSimIOModelPart)
{
    KRATOS_TRY

    KRATOS_EXPECT_TRUE(rKratosModelPart.IsDistributed());
    KRATOS_EXPECT_TRUE(rKratosModelPart.HasNodalSolutionStepVariable(PARTITION_INDEX));

    // serial checks
    CheckModelPartsAreEqual(rKratosModelPart, rCoSimIOModelPart);

    // check local nodes
    KRATOS_EXPECT_EQ(rKratosModelPart.GetCommunicator().LocalMesh().NumberOfNodes(), rCoSimIOModelPart.NumberOfLocalNodes());
    KRATOS_EXPECT_EQ(rKratosModelPart.GetCommunicator().GhostMesh().NumberOfNodes(), rCoSimIOModelPart.NumberOfGhostNodes());

    for (std::size_t i=0; i<rKratosModelPart.GetCommunicator().LocalMesh().NumberOfNodes(); ++i) {
        CheckNodesAreEqual(*(rKratosModelPart.GetCommunicator().LocalMesh().NodesBegin()+i), **(rCoSimIOModelPart.GetLocalModelPart().NodesBegin()+i));
    }

    // check ghost nodes
    // first check if the number of ghost nodes is equal
    CheckNumberOfGhostNodes(rKratosModelPart, rCoSimIOModelPart);

    // then check if the nodes are equal
    std::vector<std::size_t> ghost_nodes_per_rank_counter(rKratosModelPart.GetCommunicator().TotalProcesses(), 0);
    for (const auto& r_kratos_node : rKratosModelPart.GetCommunicator().GhostMesh().Nodes()) {
        const int part_index = r_kratos_node.FastGetSolutionStepValue(PARTITION_INDEX);
        CoSimIO::ModelPart& r_cosimio_part_mp = *rCoSimIOModelPart.GetPartitionModelParts().at(part_index);
        std::size_t rank_node_counter = ghost_nodes_per_rank_counter[part_index];
        const CoSimIO::Node& r_cosimio_node = **(r_cosimio_part_mp.NodesBegin()+rank_node_counter);

        CheckNodesAreEqual(r_kratos_node, r_cosimio_node);
        ghost_nodes_per_rank_counter[part_index]++;
    }

    KRATOS_CATCH("")
}

void CheckDistributedModelPartsAreEqualButEntitiesAreOrderedDifferently(
    const Kratos::ModelPart& rKratosModelPart,
    const CoSimIO::ModelPart& rCoSimIOModelPart)
{
    KRATOS_TRY

    KRATOS_EXPECT_TRUE(rKratosModelPart.IsDistributed());
    KRATOS_EXPECT_TRUE(rKratosModelPart.HasNodalSolutionStepVariable(PARTITION_INDEX));
    KRATOS_EXPECT_TRUE(rKratosModelPart.Has(NODES_ID_INDEX_MAP));
    KRATOS_EXPECT_TRUE(rKratosModelPart.Has(ELEMENTS_ID_INDEX_MAP));

    // serial checks
    CheckModelPartsAreEqualButEntitiesAreOrderedDifferently(rKratosModelPart, rCoSimIOModelPart);

    // check local nodes
    KRATOS_EXPECT_EQ(rKratosModelPart.GetCommunicator().LocalMesh().NumberOfNodes(), rCoSimIOModelPart.NumberOfLocalNodes());
    KRATOS_EXPECT_EQ(rKratosModelPart.GetCommunicator().GhostMesh().NumberOfNodes(), rCoSimIOModelPart.NumberOfGhostNodes());

    // check local nodes
    CheckUnorderedEntitiesAreEqual(rKratosModelPart.GetCommunicator().LocalMesh().Nodes(), rCoSimIOModelPart.LocalNodes());

    // check ghost nodes
    // first check if the number of ghost nodes is equal
    CheckNumberOfGhostNodes(rKratosModelPart, rCoSimIOModelPart);

    // then check if the nodes are equal
    CheckUnorderedEntitiesAreEqual(rKratosModelPart.GetCommunicator().GhostMesh().Nodes(), rCoSimIOModelPart.GhostNodes());

    const auto& r_kratos_ghost_nodes = rKratosModelPart.GetCommunicator().GhostMesh().Nodes();

    for (const auto& r_partition_pair : rCoSimIOModelPart.GetPartitionModelParts()) {
        const int partition_index = r_partition_pair.first;
        const auto& rp_partition_model_part = r_partition_pair.second;

        for (const auto& r_cosimio_node : rp_partition_model_part->Nodes()) {
            const auto it_kratos_node = r_kratos_ghost_nodes.find(r_cosimio_node.Id());
            KRATOS_ERROR_IF(it_kratos_node == r_kratos_ghost_nodes.end()) << "Ghost node with Id " << r_cosimio_node.Id() << " is not in the Kratos ModelPart!" << std::endl;
            const auto& r_kratos_node = *it_kratos_node;
            CheckNodesAreEqual(r_kratos_node, r_cosimio_node);
            KRATOS_EXPECT_EQ(partition_index, r_kratos_node.FastGetSolutionStepValue(PARTITION_INDEX));
        }
    }

    KRATOS_CATCH("")
}

} // namespace Kratos
} // namespace Testing
