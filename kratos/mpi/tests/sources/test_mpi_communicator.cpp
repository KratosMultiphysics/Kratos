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

#include "../applications/TrilinosApplication/custom_utilities/parallel_fill_communicator.h"

#include "testing/testing.h"

namespace Kratos {

namespace Testing {

namespace Internals {

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

KRATOS_TEST_CASE_IN_SUITE(MPICommunicatorSynchronizeOr, KratosMPICoreFastSuite)
{
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("TestModelPart");
    MPIDataCommunicator world_comm(MPI_COMM_WORLD);
    Internals::ModelPartForMPICommunicatorTests(r_model_part, world_comm);

    const int rank = world_comm.Rank();
    const int size = world_comm.Size();
    Node<3>& r_center = r_model_part.Nodes()[1];

    // Single flag
    r_center.Set(STRUCTURE, (rank == size-1));
    r_model_part.GetCommunicator().SynchronizeOrNodalFlags(STRUCTURE);
    KRATOS_CHECK_EQUAL(r_center.Is(STRUCTURE), true);

    // Multiple flags
    const bool rank_is_even( (rank % 2) == 0 );
    r_center.Clear();
    r_center.Set(INLET, rank_is_even);
    r_center.Set(OUTLET, rank_is_even);
    r_center.Set(PERIODIC, rank_is_even);

    r_model_part.GetCommunicator().SynchronizeOrNodalFlags( INLET | OUTLET );
    if (size > 1) {
        KRATOS_CHECK_EQUAL(r_center.Is(INLET), true);
        KRATOS_CHECK_EQUAL(r_center.Is(OUTLET), true);
    }
    else {
        // if there is only one rank, no communication happens
        KRATOS_CHECK_EQUAL(r_center.Is(INLET), rank_is_even);
        KRATOS_CHECK_EQUAL(r_center.Is(OUTLET), rank_is_even);
    }
    KRATOS_CHECK_EQUAL(r_center.Is(PERIODIC), rank_is_even); // This one should be left untouched

}

KRATOS_TEST_CASE_IN_SUITE(MPICommunicatorReduceAnd, KratosMPICoreFastSuite)
{
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("TestModelPart");
    MPIDataCommunicator world_comm(MPI_COMM_WORLD);
    Internals::ModelPartForMPICommunicatorTests(r_model_part, world_comm);

    const int rank = world_comm.Rank();
    const int size = world_comm.Size();
    Node<3>& r_center = r_model_part.Nodes()[1];

    // Single flag
    r_center.Set(STRUCTURE, (rank == size-1));
    r_model_part.GetCommunicator().SynchronizeAndNodalFlags(STRUCTURE);
    if (size > 1) {
        KRATOS_CHECK_EQUAL(r_center.Is(STRUCTURE), false);
    }

    // Multiple flags
    const bool rank_is_even( (rank % 2) == 0 );
    r_center.Clear();
    r_center.Set(INLET, rank_is_even);
    r_center.Set(OUTLET, rank_is_even);
    r_center.Set(PERIODIC, rank_is_even);

    r_model_part.GetCommunicator().SynchronizeAndNodalFlags( INLET | OUTLET );
    if (size > 1) {
        KRATOS_CHECK_EQUAL(r_center.Is(INLET), false);
        KRATOS_CHECK_EQUAL(r_center.Is(OUTLET), false);
    }
    else {
        // if there is only one rank, no communication happens
        KRATOS_CHECK_EQUAL(r_center.Is(INLET), rank_is_even);
        KRATOS_CHECK_EQUAL(r_center.Is(OUTLET), rank_is_even);
    }
    KRATOS_CHECK_EQUAL(r_center.Is(PERIODIC), rank_is_even); // This one should be left untouched
}

}

}