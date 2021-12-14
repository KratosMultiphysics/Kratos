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
#include "co_simulation_testing_utilities.h"


namespace Kratos {
namespace Testing {

void CheckNodesAreEqual(
    const Kratos::Node<3>& rKratosNode,
    const CoSimIO::Node& rCoSimIONode)
{
    KRATOS_TRY

    KRATOS_CHECK_EQUAL(rKratosNode.Id(), static_cast<std::size_t>(rCoSimIONode.Id()));

    KRATOS_CHECK_DOUBLE_EQUAL(rKratosNode.X(),  rCoSimIONode.X());
    KRATOS_CHECK_DOUBLE_EQUAL(rKratosNode.X0(), rCoSimIONode.X());

    KRATOS_CHECK_DOUBLE_EQUAL(rKratosNode.Y(),  rCoSimIONode.Y());
    KRATOS_CHECK_DOUBLE_EQUAL(rKratosNode.Y0(), rCoSimIONode.Y());

    KRATOS_CHECK_DOUBLE_EQUAL(rKratosNode.Z(),  rCoSimIONode.Z());
    KRATOS_CHECK_DOUBLE_EQUAL(rKratosNode.Z0(), rCoSimIONode.Z());

    KRATOS_CATCH("")
}

void CheckElementsAreEqual(
    const Kratos::Element& rKratosElement,
    const CoSimIO::Element& rCoSimIOElement)
{
    KRATOS_TRY

    // basic checks
    KRATOS_CHECK_EQUAL(rKratosElement.Id(), static_cast<std::size_t>(rCoSimIOElement.Id()));
    // KRATOS_CHECK_EQUAL(rKratosElement.Type(), rCoSimIOElement.Type());
    KRATOS_CHECK_EQUAL(rKratosElement.GetGeometry().PointsNumber(), rCoSimIOElement.NumberOfNodes());

    // check nodes
    for (std::size_t i=0; i<rCoSimIOElement.NumberOfNodes(); ++i) {
        CheckNodesAreEqual(*(rKratosElement.GetGeometry().begin()+i), **(rCoSimIOElement.NodesBegin()+i));
    }

    KRATOS_CATCH("")
}

void CheckModelPartsAreEqual(
    const Kratos::ModelPart& rKratosModelPart,
    const CoSimIO::ModelPart& rCoSimIOModelPart)
{
    KRATOS_TRY

    // basic checks
    KRATOS_CHECK_EQUAL(rKratosModelPart.NumberOfNodes(),      rCoSimIOModelPart.NumberOfNodes());
    KRATOS_CHECK_EQUAL(rKratosModelPart.NumberOfElements(),   rCoSimIOModelPart.NumberOfElements());

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

void CheckDistributedModelPartsAreEqual(
    const Kratos::ModelPart& rKratosModelPart,
    const CoSimIO::ModelPart& rCoSimIOModelPart)
{
    KRATOS_TRY

    KRATOS_CHECK(rKratosModelPart.IsDistributed());
    KRATOS_CHECK(rKratosModelPart.HasNodalSolutionStepVariable(PARTITION_INDEX));

    // serial checks
    CheckModelPartsAreEqual(rKratosModelPart, rCoSimIOModelPart);

    // check local nodes
    KRATOS_CHECK_EQUAL(rKratosModelPart.GetCommunicator().LocalMesh().NumberOfNodes(), rCoSimIOModelPart.NumberOfLocalNodes());
    KRATOS_CHECK_EQUAL(rKratosModelPart.GetCommunicator().GhostMesh().NumberOfNodes(), rCoSimIOModelPart.NumberOfGhostNodes());

    for (std::size_t i=0; i<rKratosModelPart.GetCommunicator().LocalMesh().NumberOfNodes(); ++i) {
        CheckNodesAreEqual(*(rKratosModelPart.GetCommunicator().LocalMesh().NodesBegin()+i), **(rCoSimIOModelPart.GetLocalModelPart().NodesBegin()+i));
    }

    // check ghost nodes
    // first check if the number of ghost nodes is equal
    std::vector<std::size_t> kratos_ghost_nodes_per_rank(rKratosModelPart.GetCommunicator().TotalProcesses(), 0);
    for (const auto& r_kratos_node : rKratosModelPart.GetCommunicator().GhostMesh().Nodes()) {
        kratos_ghost_nodes_per_rank[r_kratos_node.FastGetSolutionStepValue(PARTITION_INDEX)]++;
    }

    for (int part_index=0; part_index<rKratosModelPart.GetCommunicator().TotalProcesses(); ++part_index) {
        const std::size_t num_kratos_nodes_this_rank = kratos_ghost_nodes_per_rank[part_index];
        if (num_kratos_nodes_this_rank > 0) {
            KRATOS_CHECK_EQUAL(rCoSimIOModelPart.GetPartitionModelParts().count(part_index), 1);
            CoSimIO::ModelPart& r_cosimio_part_mp = *rCoSimIOModelPart.GetPartitionModelParts().at(part_index);
            KRATOS_CHECK_EQUAL(num_kratos_nodes_this_rank, r_cosimio_part_mp.NumberOfNodes());
        } else {
            KRATOS_CHECK_EQUAL(rCoSimIOModelPart.GetPartitionModelParts().count(part_index), 0);
        }
    }

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

} // namespace Kratos
} // namespace Testing
