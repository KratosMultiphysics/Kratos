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

#include <cmath>

#include "containers/model.h"
#include "includes/model_part.h"
#include "mpi/includes/mpi_communicator.h"

#include "../applications/trilinos_application/custom_utilities/parallel_fill_communicator.h"

#include "testing/testing.h"

namespace Kratos {

namespace Testing {

namespace Internals {

//void ModelPartForMPICommunicatorTests(ModelPart& rModelPart, const DataCommunicator& rComm)
//{
//    rModelPart.AddNodalSolutionStepVariable(PARTITION_INDEX);
//    Properties::Pointer p_properties = rModelPart.GetProperties(0);
//
//    const int rank = rComm.Rank();
//    const int size = rComm.Size();
//    constexpr double side_length = 1.0;
//
//    /*
//            1---3---5 ... ( i )---(i+2)
//            | \ | \ |       |   \   |
//            0---2---4 ... (i-1)---(i+1)
//    rank:    <0> <1>  ...      <n>
//    */
//    if (rank == 0)
//    {
//        rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
//    }
//
//    unsigned int index = 2*(rank + 1);
//    double x_offset = side_length*rank;
//
//    rModelPart.CreateNewNode(index  , x_offset              , side_length, 0.0);
//    rModelPart.CreateNewNode(index+1, x_offset + side_length, 0.0        , 0.0);
//
//    if (rank = size -1)
//    {
//        rModelPart.CreateNewNode(index+2, x_offset + side_length, side_length, 0.0);
//    }
//}

void ModelPartForMPICommunicatorTests(ModelPart& rModelPart, const DataCommunicator& rComm)
{
    constexpr double pi = 3.141592653589793238462643383279502884197169399375105820974944592308;
    constexpr double total_angle = pi/2.0;
    constexpr double side_length = 1.0;

    rModelPart.AddNodalSolutionStepVariable(PARTITION_INDEX);
    Properties::Pointer p_properties = rModelPart.pGetProperties(0);

    const int rank = rComm.Rank();
    const int size = rComm.Size();

    auto p_center = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    p_center->FastGetSolutionStepValue(PARTITION_INDEX) = 0.0;

    const double angle_start = rank   * (total_angle / size);
    const double angle_end   = rank+1 * (total_angle / size);

    const double x1 = side_length * std::cos(angle_start);
    const double y1 = side_length * std::sin(angle_start);

    const double x2 = side_length * std::cos(angle_end);
    const double y2 = side_length * std::sin(angle_end);

    const unsigned int index = rank + 2;
    auto p_node_1 = rModelPart.CreateNewNode(index  , x1, y1, 0.0);
    auto p_node_2 = rModelPart.CreateNewNode(index+1, x2, y2, 0.0);

    p_node_1->FastGetSolutionStepValue(PARTITION_INDEX) = 1.0*rank;
    const int remote_index = (rank == size-1) ? rank : rank + 1;
    p_node_2->FastGetSolutionStepValue(PARTITION_INDEX) = 1.0*remote_index;

    std::vector<ModelPart::IndexType> element_nodes{1, index, index+1};
    rModelPart.CreateNewElement("Element2D3N", rank+1, element_nodes, p_properties);

    ParallelFillCommunicator(rModelPart).Execute();
}

}

KRATOS_TEST_CASE_IN_SUITE(MPICommunicatorReduceOr, KratosMPICoreFastSuite)
{
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("TestModelPart");
    MPIDataCommunicator world_comm(MPI_COMM_WORLD);
    Internals::ModelPartForMPICommunicatorTests(r_model_part, world_comm);

    Node<3>& r_center = r_model_part.Nodes()[1];
    r_center.Set(STRUCTURE, (world_comm.Rank() == world_comm.Size()-1));

    r_model_part.GetCommunicator().ReduceOrNodalFlags(STRUCTURE);
    //std::cout << "Rank: " << world_comm.Rank() << ": " << r_center.Is(STRUCTURE) << "." << std::endl;
    KRATOS_CHECK_EQUAL(r_center.Is(STRUCTURE), true);
}

KRATOS_TEST_CASE_IN_SUITE(MPICommunicatorReduceAnd, KratosMPICoreFastSuite)
{
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("TestModelPart");
    MPIDataCommunicator world_comm(MPI_COMM_WORLD);
    Internals::ModelPartForMPICommunicatorTests(r_model_part, world_comm);

    Node<3>& r_center = r_model_part.Nodes()[1];
    r_center.Set(STRUCTURE, (world_comm.Rank() == world_comm.Size()-1));

    r_model_part.GetCommunicator().ReduceAndNodalFlags(STRUCTURE);
    //std::cout << "Rank: " << world_comm.Rank() << ": " << r_center.Is(STRUCTURE) << "." << std::endl;
    KRATOS_CHECK_EQUAL(r_center.Is(STRUCTURE), false);
}

}

}