//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//
//

#include "mpi.h"

#include "containers/model.h"
#include "includes/data_communicator.h"
#include "includes/model_part.h"

#include "testing/testing.h"

namespace Kratos {

namespace Testing {

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorSerializedSendRecv, KratosMPICoreFastSuite)
{
    const DataCommunicator& r_comm = DataCommunicator::GetDefault();
    const int world_rank = r_comm.Rank();
    const int world_size = r_comm.Size();
    const int send_rank = world_rank + 1 == world_size ? 0 : world_rank + 1;
    const int recv_rank = world_rank == 0 ? world_size - 1 : world_rank - 1;

    Model model;
    ModelPart& model_part = model.CreateModelPart("Send");
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.CreateNewNode(world_rank, 0.0, 0.0, 0.1*world_rank);
    for (ModelPart::NodeIterator it_node = model_part.NodesBegin(); it_node != model_part.NodesEnd(); ++it_node)
    {
        it_node->FastGetSolutionStepValue(TEMPERATURE) = 10.0*world_rank;
    }

    auto recv_node_container = r_comm.SendRecv(model_part.Nodes(), recv_rank, send_rank);

    for (auto& node: recv_node_container)
    {
        KRATOS_CHECK_EQUAL(node.Id(), (unsigned int)send_rank);
        KRATOS_CHECK_EQUAL(node.Z(), 0.1*send_rank);
        KRATOS_CHECK_EQUAL(node.FastGetSolutionStepValue(TEMPERATURE), 10.0*send_rank);
    }
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorSerializedSendAndRecv, KratosMPICoreFastSuite)
{
    const DataCommunicator& r_comm = DataCommunicator::GetDefault();
    const int world_rank = r_comm.Rank();
    const int world_size = r_comm.Size();
    const int send_rank = world_size -1;
    const int recv_rank = 0;

    if (world_rank == send_rank)
    {
        Model model;
        ModelPart& model_part = model.CreateModelPart("Send");
        model_part.AddNodalSolutionStepVariable(TEMPERATURE);
        model_part.CreateNewNode(world_rank, 0.0, 0.0, 0.1*world_rank);
        for (ModelPart::NodeIterator it_node = model_part.NodesBegin(); it_node != model_part.NodesEnd(); ++it_node)
        {
            it_node->FastGetSolutionStepValue(TEMPERATURE) = 10.0*world_rank;
        }

        r_comm.Send(model_part.Nodes(), recv_rank);
    }
    if (world_rank == recv_rank)
    {
        Model model;
        ModelPart& model_part = model.CreateModelPart("Recv");
        model_part.AddNodalSolutionStepVariable(TEMPERATURE);

        r_comm.Recv(model_part.Nodes(), send_rank);

        for (auto& node: model_part.Nodes())
        {
            KRATOS_CHECK_EQUAL(node.Id(), (unsigned int)send_rank);
            KRATOS_CHECK_EQUAL(node.Z(), 0.1*send_rank);
            KRATOS_CHECK_EQUAL(node.FastGetSolutionStepValue(TEMPERATURE), 10.0*send_rank);
        }
    }
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorSerializedBroadcast, KratosMPICoreFastSuite)
{
    const DataCommunicator& r_comm = DataCommunicator::GetDefault();
    const int world_rank = r_comm.Rank();
    const int world_size = r_comm.Size();
    const int source_rank = world_size -1;

    Model model;
    ModelPart& model_part = model.CreateModelPart("Broadcast");
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.CreateNewNode(world_rank, 0.0, 0.0, 0.1*world_rank);

    if (world_rank == source_rank)
    {
        for (ModelPart::NodeIterator it_node = model_part.NodesBegin(); it_node != model_part.NodesEnd(); ++it_node)
        {
            it_node->FastGetSolutionStepValue(TEMPERATURE) = 10.0*world_rank;
        }
    }

    r_comm.Broadcast(model_part.Nodes(), source_rank);

    for (auto& node: model_part.Nodes())
    {
        KRATOS_CHECK_EQUAL(node.Id(), (unsigned int)source_rank);
        KRATOS_CHECK_EQUAL(node.Z(), 0.1*source_rank);
        KRATOS_CHECK_EQUAL(node.FastGetSolutionStepValue(TEMPERATURE), 10.0*source_rank);
    }
}


}

}