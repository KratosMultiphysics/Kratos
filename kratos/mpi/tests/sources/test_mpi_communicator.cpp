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
#include "includes/parallel_environment.h"
#include "mpi/includes/mpi_communicator.h"
#include "mpi/utilities/parallel_fill_communicator.h"

#include "testing/testing.h"

namespace Kratos {

namespace Testing {

namespace Internals {

void ModelPartForMPICommunicatorTests(ModelPart& rModelPart, const DataCommunicator& rComm)
{
    /* NOTE: the modelpart should at least have PARTITION_INDEX in the nodal solution step data */
    constexpr double pi = 3.141592653589793238462643383279502884197169399375105820974944592308;
    constexpr double total_angle = pi/2.0;
    constexpr double side_length = 1.0;

    Properties::Pointer p_properties = rModelPart.CreateNewProperties(0);

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

    const unsigned int local_index = rank + 2;
    const unsigned int ghost_index = (size == 1) || (rank != size-1) ? rank + 3 : 2;
    auto p_node_1 = rModelPart.CreateNewNode(local_index, x1, y1, 0.0);
    auto p_node_2 = rModelPart.CreateNewNode(ghost_index, x2, y2, 0.0);

    p_node_1->FastGetSolutionStepValue(PARTITION_INDEX) = 1.0*rank;
    const int remote_rank = (rank != size-1) ? rank + 1 : 0;
    p_node_2->FastGetSolutionStepValue(PARTITION_INDEX) = 1.0*remote_rank;

    std::vector<ModelPart::IndexType> element_nodes{1, local_index, ghost_index};
    rModelPart.CreateNewElement("Element2D3N", rank+1, element_nodes, p_properties);

    ParallelFillCommunicator(rModelPart).Execute();
}

}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPICommunicatorCreation, KratosMPICoreFastSuite)
{
    const DataCommunicator& r_world = ParallelEnvironment::GetDataCommunicator("World");

    Model model;
    ModelPart& r_model_part = model.CreateModelPart("TestModelPart");
    r_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    // Test regular construction
    const MPICommunicator constructed_communicator = MPICommunicator(&r_model_part.GetNodalSolutionStepVariablesList(), r_world);

    KRATOS_CHECK_EQUAL(constructed_communicator.IsDistributed(), true);
    KRATOS_CHECK_EQUAL(constructed_communicator.MyPID(), r_world.Rank());
    KRATOS_CHECK_EQUAL(constructed_communicator.TotalProcesses(), r_world.Size());

    // Test creation with given DataCommunicator
    Communicator::Pointer p_created_communicator = constructed_communicator.Create(r_world);

    KRATOS_CHECK_EQUAL(p_created_communicator->IsDistributed(), true);
    KRATOS_CHECK_EQUAL(p_created_communicator->MyPID(), r_world.Rank());
    KRATOS_CHECK_EQUAL(p_created_communicator->TotalProcesses(), r_world.Size());

    // Test creation using reference's DataCommunicator
    p_created_communicator.reset();
    p_created_communicator = constructed_communicator.Create();

    KRATOS_CHECK_EQUAL(p_created_communicator->IsDistributed(), true);
    KRATOS_CHECK_EQUAL(p_created_communicator->MyPID(), r_world.Rank());
    KRATOS_CHECK_EQUAL(p_created_communicator->TotalProcesses(), r_world.Size());
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPICommunicatorSynchronizeOr, KratosMPICoreFastSuite)
{
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("TestModelPart");
    r_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

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

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPICommunicatorSynchronizeAnd, KratosMPICoreFastSuite)
{
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("TestModelPart");
    r_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

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

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPICommunicatorSynchronizeNodalFlags, KratosMPICoreFastSuite)
{
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("TestModelPart");
    r_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    MPIDataCommunicator world_comm(MPI_COMM_WORLD);
    Internals::ModelPartForMPICommunicatorTests(r_model_part, world_comm);

    const int rank = world_comm.Rank();
    const bool rank_is_even( (rank % 2) == 0 );

    for (auto i_node = r_model_part.NodesBegin(); i_node != r_model_part.NodesEnd(); ++i_node)
    {
        i_node->Set(INLET, rank_is_even);
        i_node->Set(OUTLET, rank_is_even);
        i_node->Set(PERIODIC, !rank_is_even);
    }

    r_model_part.GetCommunicator().SynchronizeNodalFlags();
    // End result: the entire Flags are copied from owner rank to ghost copies (both defined status and values)
    for (auto i_node = r_model_part.NodesBegin(); i_node != r_model_part.NodesEnd(); ++i_node)
    {
        int owner_rank = i_node->FastGetSolutionStepValue(PARTITION_INDEX, 0);
        bool owner_is_even = ((owner_rank % 2) == 0);
        KRATOS_CHECK_EQUAL(i_node->Is(INLET), owner_is_even);
        KRATOS_CHECK_EQUAL(i_node->Is(OUTLET), owner_is_even);
        KRATOS_CHECK_EQUAL(i_node->Is(PERIODIC), !owner_is_even);
        KRATOS_CHECK_EQUAL(i_node->IsDefined(SLIP), false); // this one was never set
    }
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPICommunicatorNodalSolutionStepVariableAssembly, KratosMPICoreFastSuite)
{
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("TestModelPart");
    r_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);
    r_model_part.AddNodalSolutionStepVariable(DOMAIN_SIZE);           // Variable<int>
    r_model_part.AddNodalSolutionStepVariable(TEMPERATURE);           // Variable<double>
    r_model_part.AddNodalSolutionStepVariable(VELOCITY);              // Variable< array_1d<double,3> >
    r_model_part.AddNodalSolutionStepVariable(CAUCHY_STRESS_VECTOR);  // Variable<Vector>
    r_model_part.AddNodalSolutionStepVariable(DEFORMATION_GRADIENT);  // Variable<Matrix>

    MPIDataCommunicator comm_world(MPI_COMM_WORLD);
    Internals::ModelPartForMPICommunicatorTests(r_model_part, comm_world);
    int rank = comm_world.Rank();
    int size = comm_world.Size();

    for (auto i_node = r_model_part.NodesBegin(); i_node != r_model_part.NodesEnd(); ++i_node)
    {
        i_node->FastGetSolutionStepValue(DOMAIN_SIZE, 0) = 1;
        i_node->FastGetSolutionStepValue(TEMPERATURE, 0) = 2.0;
        i_node->FastGetSolutionStepValue(VELOCITY_X, 0) = 1.0;
        i_node->FastGetSolutionStepValue(VELOCITY_Y, 0) = 2.0;
        Vector& r_cauchy_stress = i_node->FastGetSolutionStepValue(CAUCHY_STRESS_VECTOR, 0);
        r_cauchy_stress.resize(2, false);
        r_cauchy_stress = ZeroVector(2);
        r_cauchy_stress[1] = 1.0;
        Matrix& r_deformation_gradient = i_node->FastGetSolutionStepValue(DEFORMATION_GRADIENT, 0);
        r_deformation_gradient.resize(3,2,false);
        r_deformation_gradient = ZeroMatrix(3,2);
        r_deformation_gradient(2,0) = 1.0;
    }

    Communicator& r_comm = r_model_part.GetCommunicator();

    // center is local to rank 0 and ghost in all other ranks
    Node<3>& r_center = r_model_part.Nodes()[1];
    // local and ghost nodes are each known in two ranks
    const unsigned int local_id = rank + 2;
    unsigned int ghost_id = (size == 1) || (rank != size-1) ? rank + 3 : 2;
    Node<3>& r_local = r_model_part.Nodes()[local_id];
    Node<3>& r_ghost = r_model_part.Nodes()[ghost_id];

    r_comm.AssembleCurrentData(DOMAIN_SIZE);
    int expected_int = (size > 1) ? 2 : 1;
    KRATOS_CHECK_EQUAL(r_center.FastGetSolutionStepValue(DOMAIN_SIZE,0), size);
    KRATOS_CHECK_EQUAL( r_local.FastGetSolutionStepValue(DOMAIN_SIZE,0), expected_int);
    KRATOS_CHECK_EQUAL( r_ghost.FastGetSolutionStepValue(DOMAIN_SIZE,0), expected_int);

    r_comm.AssembleCurrentData(TEMPERATURE);
    int expected_double = 2.0*expected_int;
    KRATOS_CHECK_EQUAL(r_center.FastGetSolutionStepValue(TEMPERATURE,0), 2.0*size);
    KRATOS_CHECK_EQUAL( r_local.FastGetSolutionStepValue(TEMPERATURE,0), expected_double);
    KRATOS_CHECK_EQUAL( r_ghost.FastGetSolutionStepValue(TEMPERATURE,0), expected_double);

    r_comm.AssembleCurrentData(VELOCITY);
    KRATOS_CHECK_EQUAL(r_center.FastGetSolutionStepValue(VELOCITY_X,0), 1.0*size);
    KRATOS_CHECK_EQUAL( r_local.FastGetSolutionStepValue(VELOCITY_X,0), 1.0*expected_int);
    KRATOS_CHECK_EQUAL( r_ghost.FastGetSolutionStepValue(VELOCITY_X,0), 1.0*expected_int);

    KRATOS_CHECK_EQUAL(r_center.FastGetSolutionStepValue(VELOCITY_Y,0), 2.0*size);
    KRATOS_CHECK_EQUAL( r_local.FastGetSolutionStepValue(VELOCITY_Y,0), 2.0*expected_int);
    KRATOS_CHECK_EQUAL( r_ghost.FastGetSolutionStepValue(VELOCITY_Y,0), 2.0*expected_int);

    KRATOS_CHECK_EQUAL(r_center.FastGetSolutionStepValue(VELOCITY_Z,0), 0.0);
    KRATOS_CHECK_EQUAL( r_local.FastGetSolutionStepValue(VELOCITY_Z,0), 0.0);
    KRATOS_CHECK_EQUAL( r_ghost.FastGetSolutionStepValue(VELOCITY_Z,0), 0.0);

    r_comm.AssembleCurrentData(CAUCHY_STRESS_VECTOR);
    const auto& r_assembled_center_vector = r_center.FastGetSolutionStepValue(CAUCHY_STRESS_VECTOR,0);
    KRATOS_CHECK_EQUAL(r_assembled_center_vector.size(), 2);
    KRATOS_CHECK_EQUAL(r_assembled_center_vector[0], 0.0);
    KRATOS_CHECK_EQUAL(r_assembled_center_vector[1], 1.0*size);

    const auto& r_assembled_local_vector = r_local.FastGetSolutionStepValue(CAUCHY_STRESS_VECTOR,0);
    KRATOS_CHECK_EQUAL(r_assembled_local_vector.size(), 2);
    KRATOS_CHECK_EQUAL(r_assembled_local_vector[0], 0.0);
    KRATOS_CHECK_EQUAL(r_assembled_local_vector[1], 1.0*expected_int);

    const auto& r_assembled_ghost_vector = r_ghost.FastGetSolutionStepValue(CAUCHY_STRESS_VECTOR,0);
    KRATOS_CHECK_EQUAL(r_assembled_ghost_vector.size(), 2);
    KRATOS_CHECK_EQUAL(r_assembled_ghost_vector[0], 0.0);
    KRATOS_CHECK_EQUAL(r_assembled_ghost_vector[1], 1.0*expected_int);

    r_comm.AssembleCurrentData(DEFORMATION_GRADIENT);
    const auto& r_assembled_center_matrix = r_center.FastGetSolutionStepValue(DEFORMATION_GRADIENT,0);
    KRATOS_CHECK_EQUAL(r_assembled_center_matrix.size1(), 3);
    KRATOS_CHECK_EQUAL(r_assembled_center_matrix.size2(), 2);
    KRATOS_CHECK_EQUAL(r_assembled_center_matrix(0,0), 0.0);
    KRATOS_CHECK_EQUAL(r_assembled_center_matrix(2,0), 1.0*size);

    const auto& r_assembled_local_matrix = r_local.FastGetSolutionStepValue(DEFORMATION_GRADIENT,0);
    KRATOS_CHECK_EQUAL(r_assembled_local_matrix.size1(), 3);
    KRATOS_CHECK_EQUAL(r_assembled_local_matrix.size2(), 2);
    KRATOS_CHECK_EQUAL(r_assembled_local_matrix(0,0), 0.0);
    KRATOS_CHECK_EQUAL(r_assembled_local_matrix(2,0), 1.0*expected_int);

    const auto& r_assembled_ghost_matrix = r_ghost.FastGetSolutionStepValue(DEFORMATION_GRADIENT,0);
    KRATOS_CHECK_EQUAL(r_assembled_ghost_matrix.size1(), 3);
    KRATOS_CHECK_EQUAL(r_assembled_ghost_matrix.size2(), 2);
    KRATOS_CHECK_EQUAL(r_assembled_ghost_matrix(0,0), 0.0);
    KRATOS_CHECK_EQUAL(r_assembled_ghost_matrix(2,0), 1.0*expected_int);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPICommunicatorNodalSolutionStepVariableSynchronize, KratosMPICoreFastSuite)
{
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("TestModelPart");
    r_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);
    r_model_part.AddNodalSolutionStepVariable(DOMAIN_SIZE);           // Variable<int>
    r_model_part.AddNodalSolutionStepVariable(TEMPERATURE);           // Variable<double>
    r_model_part.AddNodalSolutionStepVariable(IS_RESTARTED);          // Variable<bool>
    r_model_part.AddNodalSolutionStepVariable(VELOCITY);              // Variable< array_1d<double,3> >
    r_model_part.AddNodalSolutionStepVariable(CAUCHY_STRESS_VECTOR);  // Variable<Vector>
    r_model_part.AddNodalSolutionStepVariable(DEFORMATION_GRADIENT);  // Variable<Matrix>

    MPIDataCommunicator comm_world(MPI_COMM_WORLD);
    Internals::ModelPartForMPICommunicatorTests(r_model_part, comm_world);

    Communicator& r_comm = r_model_part.GetCommunicator();
    ModelPart::NodesContainerType& r_local_nodes = r_comm.LocalMesh().Nodes();
    for (auto i_node = r_local_nodes.begin(); i_node != r_local_nodes.end(); ++i_node)
    {
        i_node->FastGetSolutionStepValue(DOMAIN_SIZE, 0) = 1;
        i_node->FastGetSolutionStepValue(TEMPERATURE, 0) = 2.0;
        i_node->FastGetSolutionStepValue(IS_RESTARTED, 0) = true;
        i_node->FastGetSolutionStepValue(VELOCITY_X, 0) = 1.0;
        i_node->FastGetSolutionStepValue(VELOCITY_Y, 0) = 2.0;
        Vector& r_cauchy_stress = i_node->FastGetSolutionStepValue(CAUCHY_STRESS_VECTOR, 0);
        r_cauchy_stress.resize(2, false);
        r_cauchy_stress = ZeroVector(2);
        r_cauchy_stress[1] = 1.0;
        Matrix& r_deformation_gradient = i_node->FastGetSolutionStepValue(DEFORMATION_GRADIENT, 0);
        r_deformation_gradient.resize(3,2,false);
        r_deformation_gradient = ZeroMatrix(3,2);
        r_deformation_gradient(2,1) = 1.0;
    }

    r_comm.SynchronizeVariable(DOMAIN_SIZE);
    for (auto i_node = r_model_part.NodesBegin(); i_node != r_model_part.NodesEnd(); ++i_node)
    {
        KRATOS_CHECK_EQUAL(i_node->FastGetSolutionStepValue(DOMAIN_SIZE,0), 1);
    }

    r_comm.SynchronizeVariable(TEMPERATURE);
    for (auto i_node = r_model_part.NodesBegin(); i_node != r_model_part.NodesEnd(); ++i_node)
    {
        KRATOS_CHECK_EQUAL(i_node->FastGetSolutionStepValue(TEMPERATURE,0), 2.0);
    }

    r_comm.SynchronizeVariable(IS_RESTARTED);
    for (auto i_node = r_model_part.NodesBegin(); i_node != r_model_part.NodesEnd(); ++i_node)
    {
        KRATOS_CHECK_EQUAL(i_node->FastGetSolutionStepValue(IS_RESTARTED,0), true);
    }

    r_comm.SynchronizeVariable(VELOCITY);
    for (auto i_node = r_model_part.NodesBegin(); i_node != r_model_part.NodesEnd(); ++i_node)
    {
        KRATOS_CHECK_EQUAL(i_node->FastGetSolutionStepValue(VELOCITY_X,0), 1.0);
        KRATOS_CHECK_EQUAL(i_node->FastGetSolutionStepValue(VELOCITY_Y,0), 2.0);
        KRATOS_CHECK_EQUAL(i_node->FastGetSolutionStepValue(VELOCITY_Z,0), 0.0);
    }

    r_comm.SynchronizeVariable(CAUCHY_STRESS_VECTOR);
    for (auto i_node = r_model_part.NodesBegin(); i_node != r_model_part.NodesEnd(); ++i_node)
    {
        const auto& r_vector = i_node->FastGetSolutionStepValue(CAUCHY_STRESS_VECTOR,0);
        KRATOS_CHECK_EQUAL(r_vector.size(), 2);
        KRATOS_CHECK_EQUAL(r_vector[0], 0.0);
        KRATOS_CHECK_EQUAL(r_vector[1], 1.0);
    }

    r_comm.SynchronizeVariable(DEFORMATION_GRADIENT);
    for (auto i_node = r_model_part.NodesBegin(); i_node != r_model_part.NodesEnd(); ++i_node)
    {
        const auto& r_matrix = i_node->FastGetSolutionStepValue(DEFORMATION_GRADIENT,0);
        KRATOS_CHECK_EQUAL(r_matrix.size1(), 3);
        KRATOS_CHECK_EQUAL(r_matrix.size2(), 2);
        KRATOS_CHECK_EQUAL(r_matrix(0,0), 0.0);
        KRATOS_CHECK_EQUAL(r_matrix(2,1), 1.0);
    }
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPICommunicatorNodalDataAssembly, KratosMPICoreFastSuite)
{
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("TestModelPart");
    r_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    MPIDataCommunicator comm_world(MPI_COMM_WORLD);
    Internals::ModelPartForMPICommunicatorTests(r_model_part, comm_world);
    int rank = comm_world.Rank();
    int size = comm_world.Size();

    for (auto i_node = r_model_part.NodesBegin(); i_node != r_model_part.NodesEnd(); ++i_node)
    {
        i_node->GetValue(DOMAIN_SIZE) = 1;
        i_node->GetValue(TEMPERATURE) = 2.0;
        i_node->GetValue(VELOCITY_X) = 1.0;
        i_node->GetValue(VELOCITY_Y) = 2.0;
        Vector& r_cauchy_stress = i_node->GetValue(CAUCHY_STRESS_VECTOR);
        r_cauchy_stress.resize(2, false);
        r_cauchy_stress = ZeroVector(2);
        r_cauchy_stress[1] = 1.0;
        Matrix& r_deformation_gradient = i_node->GetValue(DEFORMATION_GRADIENT);
        r_deformation_gradient.resize(3,2,false);
        r_deformation_gradient = ZeroMatrix(3,2);
        r_deformation_gradient(2,0) = 1.0;
    }

    Communicator& r_comm = r_model_part.GetCommunicator();

    // center is local to rank 0 and ghost in all other ranks
    Node<3>& r_center = r_model_part.Nodes()[1];
    // local and ghost nodes are each known in two ranks
    const unsigned int local_id = rank + 2;
    unsigned int ghost_id = (size == 1) || (rank != size-1) ? rank + 3 : 2;
    Node<3>& r_local = r_model_part.Nodes()[local_id];
    Node<3>& r_ghost = r_model_part.Nodes()[ghost_id];

    r_comm.AssembleNonHistoricalData(DOMAIN_SIZE);
    int expected_int = (size > 1) ? 2 : 1;
    KRATOS_CHECK_EQUAL(r_center.GetValue(DOMAIN_SIZE), size);
    KRATOS_CHECK_EQUAL( r_local.GetValue(DOMAIN_SIZE), expected_int);
    KRATOS_CHECK_EQUAL( r_ghost.GetValue(DOMAIN_SIZE), expected_int);

    r_comm.AssembleNonHistoricalData(TEMPERATURE);
    int expected_double = 2.0*expected_int;
    KRATOS_CHECK_EQUAL(r_center.GetValue(TEMPERATURE), 2.0*size);
    KRATOS_CHECK_EQUAL( r_local.GetValue(TEMPERATURE), expected_double);
    KRATOS_CHECK_EQUAL( r_ghost.GetValue(TEMPERATURE), expected_double);

    r_comm.AssembleNonHistoricalData(VELOCITY);
    KRATOS_CHECK_EQUAL(r_center.GetValue(VELOCITY_X), 1.0*size);
    KRATOS_CHECK_EQUAL( r_local.GetValue(VELOCITY_X), 1.0*expected_int);
    KRATOS_CHECK_EQUAL( r_ghost.GetValue(VELOCITY_X), 1.0*expected_int);

    KRATOS_CHECK_EQUAL(r_center.GetValue(VELOCITY_Y), 2.0*size);
    KRATOS_CHECK_EQUAL( r_local.GetValue(VELOCITY_Y), 2.0*expected_int);
    KRATOS_CHECK_EQUAL( r_ghost.GetValue(VELOCITY_Y), 2.0*expected_int);

    KRATOS_CHECK_EQUAL(r_center.GetValue(VELOCITY_Z), 0.0);
    KRATOS_CHECK_EQUAL( r_local.GetValue(VELOCITY_Z), 0.0);
    KRATOS_CHECK_EQUAL( r_ghost.GetValue(VELOCITY_Z), 0.0);

    r_comm.AssembleNonHistoricalData(CAUCHY_STRESS_VECTOR);
    const auto& r_assembled_center_vector = r_center.GetValue(CAUCHY_STRESS_VECTOR);
    KRATOS_CHECK_EQUAL(r_assembled_center_vector.size(), 2);
    KRATOS_CHECK_EQUAL(r_assembled_center_vector[0], 0.0);
    KRATOS_CHECK_EQUAL(r_assembled_center_vector[1], 1.0*size);

    const auto& r_assembled_local_vector = r_local.GetValue(CAUCHY_STRESS_VECTOR);
    KRATOS_CHECK_EQUAL(r_assembled_local_vector.size(), 2);
    KRATOS_CHECK_EQUAL(r_assembled_local_vector[0], 0.0);
    KRATOS_CHECK_EQUAL(r_assembled_local_vector[1], 1.0*expected_int);

    const auto& r_assembled_ghost_vector = r_ghost.GetValue(CAUCHY_STRESS_VECTOR);
    KRATOS_CHECK_EQUAL(r_assembled_ghost_vector.size(), 2);
    KRATOS_CHECK_EQUAL(r_assembled_ghost_vector[0], 0.0);
    KRATOS_CHECK_EQUAL(r_assembled_ghost_vector[1], 1.0*expected_int);

    r_comm.AssembleNonHistoricalData(DEFORMATION_GRADIENT);
    const auto& r_assembled_center_matrix = r_center.GetValue(DEFORMATION_GRADIENT);
    KRATOS_CHECK_EQUAL(r_assembled_center_matrix.size1(), 3);
    KRATOS_CHECK_EQUAL(r_assembled_center_matrix.size2(), 2);
    KRATOS_CHECK_EQUAL(r_assembled_center_matrix(0,0), 0.0);
    KRATOS_CHECK_EQUAL(r_assembled_center_matrix(2,0), 1.0*size);

    const auto& r_assembled_local_matrix = r_local.GetValue(DEFORMATION_GRADIENT);
    KRATOS_CHECK_EQUAL(r_assembled_local_matrix.size1(), 3);
    KRATOS_CHECK_EQUAL(r_assembled_local_matrix.size2(), 2);
    KRATOS_CHECK_EQUAL(r_assembled_local_matrix(0,0), 0.0);
    KRATOS_CHECK_EQUAL(r_assembled_local_matrix(2,0), 1.0*expected_int);

    const auto& r_assembled_ghost_matrix = r_ghost.GetValue(DEFORMATION_GRADIENT);
    KRATOS_CHECK_EQUAL(r_assembled_ghost_matrix.size1(), 3);
    KRATOS_CHECK_EQUAL(r_assembled_ghost_matrix.size2(), 2);
    KRATOS_CHECK_EQUAL(r_assembled_ghost_matrix(0,0), 0.0);
    KRATOS_CHECK_EQUAL(r_assembled_ghost_matrix(2,0), 1.0*expected_int);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPICommunicatorNodalDataSynchronize, KratosMPICoreFastSuite)
{
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("TestModelPart");
    r_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    MPIDataCommunicator comm_world(MPI_COMM_WORLD);
    Internals::ModelPartForMPICommunicatorTests(r_model_part, comm_world);

    Communicator& r_comm = r_model_part.GetCommunicator();
    ModelPart::NodesContainerType& r_local_nodes = r_comm.LocalMesh().Nodes();
    for (auto i_node = r_local_nodes.begin(); i_node != r_local_nodes.end(); ++i_node)
    {
        i_node->GetValue(DOMAIN_SIZE) = 1;
        i_node->GetValue(TEMPERATURE) = 2.0;
        i_node->GetValue(IS_RESTARTED) = true;
        i_node->GetValue(VELOCITY_X) = 1.0;
        i_node->GetValue(VELOCITY_Y) = 2.0;
        Vector& r_cauchy_stress = i_node->GetValue(CAUCHY_STRESS_VECTOR);
        r_cauchy_stress.resize(2, false);
        r_cauchy_stress = ZeroVector(2);
        r_cauchy_stress[1] = 1.0;
        Matrix& r_deformation_gradient = i_node->GetValue(DEFORMATION_GRADIENT);
        r_deformation_gradient.resize(3,2,false);
        r_deformation_gradient = ZeroMatrix(3,2);
        r_deformation_gradient(2,1) = 1.0;
    }

    r_comm.SynchronizeNonHistoricalVariable(DOMAIN_SIZE);
    for (auto i_node = r_model_part.NodesBegin(); i_node != r_model_part.NodesEnd(); ++i_node)
    {
        KRATOS_CHECK_EQUAL(i_node->GetValue(DOMAIN_SIZE), 1);
    }

    r_comm.SynchronizeNonHistoricalVariable(TEMPERATURE);
    for (auto i_node = r_model_part.NodesBegin(); i_node != r_model_part.NodesEnd(); ++i_node)
    {
        KRATOS_CHECK_EQUAL(i_node->GetValue(TEMPERATURE), 2.0);
    }

    r_comm.SynchronizeNonHistoricalVariable(IS_RESTARTED);
    for (auto i_node = r_model_part.NodesBegin(); i_node != r_model_part.NodesEnd(); ++i_node)
    {
        KRATOS_CHECK_EQUAL(i_node->GetValue(IS_RESTARTED), true);
    }

    r_comm.SynchronizeNonHistoricalVariable(VELOCITY);
    for (auto i_node = r_model_part.NodesBegin(); i_node != r_model_part.NodesEnd(); ++i_node)
    {
        KRATOS_CHECK_EQUAL(i_node->GetValue(VELOCITY_X), 1.0);
        KRATOS_CHECK_EQUAL(i_node->GetValue(VELOCITY_Y), 2.0);
        KRATOS_CHECK_EQUAL(i_node->GetValue(VELOCITY_Z), 0.0);
    }

    r_comm.SynchronizeNonHistoricalVariable(CAUCHY_STRESS_VECTOR);
    for (auto i_node = r_model_part.NodesBegin(); i_node != r_model_part.NodesEnd(); ++i_node)
    {
        const auto& r_vector = i_node->GetValue(CAUCHY_STRESS_VECTOR);
        KRATOS_CHECK_EQUAL(r_vector.size(), 2);
        KRATOS_CHECK_EQUAL(r_vector[0], 0.0);
        KRATOS_CHECK_EQUAL(r_vector[1], 1.0);
    }

    r_comm.SynchronizeNonHistoricalVariable(DEFORMATION_GRADIENT);
    for (auto i_node = r_model_part.NodesBegin(); i_node != r_model_part.NodesEnd(); ++i_node)
    {
        const auto& r_matrix = i_node->GetValue(DEFORMATION_GRADIENT);
        KRATOS_CHECK_EQUAL(r_matrix.size1(), 3);
        KRATOS_CHECK_EQUAL(r_matrix.size2(), 2);
        KRATOS_CHECK_EQUAL(r_matrix(0,0), 0.0);
        KRATOS_CHECK_EQUAL(r_matrix(2,1), 1.0);
    }
}


KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPICommunicatorNodalSolutionStepVariableSyncToMin, KratosMPICoreFastSuite)
{
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("TestModelPart");
    r_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);
    r_model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    MPIDataCommunicator comm_world(MPI_COMM_WORLD);
    Internals::ModelPartForMPICommunicatorTests(r_model_part, comm_world);
    int rank = comm_world.Rank();
    int size = comm_world.Size();

    for (auto i_node = r_model_part.NodesBegin(); i_node != r_model_part.NodesEnd(); ++i_node)
    {
        i_node->FastGetSolutionStepValue(TEMPERATURE, 0) = 10.0*rank;
    }

    Communicator& r_comm = r_model_part.GetCommunicator();

    // center is local to rank 0 and ghost in all other ranks
    Node<3>& r_center = r_model_part.Nodes()[1];
    // local and ghost nodes are each known in two ranks
    const unsigned int local_id = rank + 2;
    unsigned int ghost_id = (size == 1) || (rank != size-1) ? rank + 3 : 2;
    Node<3>& r_local = r_model_part.Nodes()[local_id];
    Node<3>& r_ghost = r_model_part.Nodes()[ghost_id];

    int expected_local = (rank > 0) ? 10.0*(rank-1) : 0.0;
    int expected_ghost = (rank + 1 < size) ? 10.0*rank : 0.0;

    r_comm.SynchronizeCurrentDataToMin(TEMPERATURE);
    KRATOS_CHECK_EQUAL(r_center.FastGetSolutionStepValue(TEMPERATURE,0), 0.0);
    KRATOS_CHECK_EQUAL( r_local.FastGetSolutionStepValue(TEMPERATURE,0), expected_local);
    KRATOS_CHECK_EQUAL( r_ghost.FastGetSolutionStepValue(TEMPERATURE,0), expected_ghost);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPICommunicatorNodalDataariableSyncToMin, KratosMPICoreFastSuite)
{
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("TestModelPart");
    r_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    MPIDataCommunicator comm_world(MPI_COMM_WORLD);
    Internals::ModelPartForMPICommunicatorTests(r_model_part, comm_world);
    int rank = comm_world.Rank();
    int size = comm_world.Size();

    for (auto i_node = r_model_part.NodesBegin(); i_node != r_model_part.NodesEnd(); ++i_node)
    {
        i_node->SetValue(TEMPERATURE, 10.0*rank);
    }

    Communicator& r_comm = r_model_part.GetCommunicator();

    // center is local to rank 0 and ghost in all other ranks
    Node<3>& r_center = r_model_part.Nodes()[1];
    // local and ghost nodes are each known in two ranks
    const unsigned int local_id = rank + 2;
    unsigned int ghost_id = (size == 1) || (rank != size-1) ? rank + 3 : 2;
    Node<3>& r_local = r_model_part.Nodes()[local_id];
    Node<3>& r_ghost = r_model_part.Nodes()[ghost_id];

    int expected_local = (rank > 0) ? 10.0*(rank-1) : 0.0;
    int expected_ghost = (rank + 1 < size) ? 10.0*rank : 0.0;

    r_comm.SynchronizeNonHistoricalDataToMin(TEMPERATURE);
    KRATOS_CHECK_EQUAL(r_center.GetValue(TEMPERATURE), 0.0);
    KRATOS_CHECK_EQUAL( r_local.GetValue(TEMPERATURE), expected_local);
    KRATOS_CHECK_EQUAL( r_ghost.GetValue(TEMPERATURE), expected_ghost);
}


KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPICommunicatorNodalSolutionStepDataSynchronize, KratosMPICoreFastSuite)
{
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("TestModelPart");
    r_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);
    r_model_part.AddNodalSolutionStepVariable(DOMAIN_SIZE);           // Variable<int>
    r_model_part.AddNodalSolutionStepVariable(TEMPERATURE);           // Variable<double>
    r_model_part.AddNodalSolutionStepVariable(IS_RESTARTED);          // Variable<bool>
    r_model_part.AddNodalSolutionStepVariable(VELOCITY);              // Variable< array_1d<double,3> >
    r_model_part.AddNodalSolutionStepVariable(CAUCHY_STRESS_VECTOR);  // Variable<Vector>
    r_model_part.AddNodalSolutionStepVariable(DEFORMATION_GRADIENT);  // Variable<Matrix>

    r_model_part.SetBufferSize(2);
    MPIDataCommunicator comm_world(MPI_COMM_WORLD);
    Internals::ModelPartForMPICommunicatorTests(r_model_part, comm_world);

    Communicator& r_comm = r_model_part.GetCommunicator();
    ModelPart::NodesContainerType& r_local_nodes = r_comm.LocalMesh().Nodes();
    for (auto i_node = r_local_nodes.begin(); i_node != r_local_nodes.end(); ++i_node)
    {
        i_node->FastGetSolutionStepValue(DOMAIN_SIZE, 0) = 1;
        i_node->FastGetSolutionStepValue(TEMPERATURE, 0) = 2.0;
        i_node->FastGetSolutionStepValue(IS_RESTARTED, 0) = true;
        i_node->FastGetSolutionStepValue(VELOCITY_X, 0) = 1.0;
        i_node->FastGetSolutionStepValue(VELOCITY_Y, 0) = 2.0;
        Vector& r_cauchy_stress = i_node->FastGetSolutionStepValue(CAUCHY_STRESS_VECTOR, 0);
        r_cauchy_stress.resize(2, false);
        r_cauchy_stress = ZeroVector(2);
        r_cauchy_stress[1] = 1.0;
        Matrix& r_deformation_gradient = i_node->FastGetSolutionStepValue(DEFORMATION_GRADIENT, 0);
        r_deformation_gradient.resize(3,2,false);
        r_deformation_gradient = ZeroMatrix(3,2);
        r_deformation_gradient(2,1) = 1.0;
    }

    // NodalSolutionStepData synchronization should also update old time steps.
    r_model_part.CloneTimeStep(1.0);

    r_comm.SynchronizeNodalSolutionStepsData();

    for (auto i_node = r_model_part.NodesBegin(); i_node != r_model_part.NodesEnd(); ++i_node)
    {
        for (unsigned int step = 0; step < 2; step++)
        {
            KRATOS_CHECK_EQUAL(i_node->FastGetSolutionStepValue(DOMAIN_SIZE,step), 1);
            KRATOS_CHECK_EQUAL(i_node->FastGetSolutionStepValue(TEMPERATURE,step), 2.0);
            KRATOS_CHECK_EQUAL(i_node->FastGetSolutionStepValue(IS_RESTARTED,step), true);
            KRATOS_CHECK_EQUAL(i_node->FastGetSolutionStepValue(VELOCITY_X,step), 1.0);
            KRATOS_CHECK_EQUAL(i_node->FastGetSolutionStepValue(VELOCITY_Y,step), 2.0);
            KRATOS_CHECK_EQUAL(i_node->FastGetSolutionStepValue(VELOCITY_Z,step), 0.0);
            const auto& r_vector = i_node->FastGetSolutionStepValue(CAUCHY_STRESS_VECTOR,step);
            KRATOS_CHECK_EQUAL(r_vector.size(), 2);
            KRATOS_CHECK_EQUAL(r_vector[0], 0.0);
            KRATOS_CHECK_EQUAL(r_vector[1], 1.0);
            const auto& r_matrix = i_node->FastGetSolutionStepValue(DEFORMATION_GRADIENT,step);
            KRATOS_CHECK_EQUAL(r_matrix.size1(), 3);
            KRATOS_CHECK_EQUAL(r_matrix.size2(), 2);
            KRATOS_CHECK_EQUAL(r_matrix(0,0), 0.0);
            KRATOS_CHECK_EQUAL(r_matrix(2,1), 1.0);
        }
    }
}


KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPICommunicatorSynchronizeDofIds, KratosMPICoreFastSuite)
{
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("TestModelPart");
    r_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);
    r_model_part.AddNodalSolutionStepVariable(VELOCITY);
    r_model_part.AddNodalSolutionStepVariable(REACTION);

    MPIDataCommunicator comm_world(MPI_COMM_WORLD);
    Internals::ModelPartForMPICommunicatorTests(r_model_part, comm_world);

    for (auto i_node = r_model_part.NodesBegin(); i_node != r_model_part.NodesEnd(); ++i_node)
    {
        i_node->AddDof(VELOCITY_X, REACTION_X);
        i_node->AddDof(VELOCITY_Y, REACTION_Y);
        i_node->AddDof(VELOCITY_Z, REACTION_Z);
    }

    Communicator& r_comm = r_model_part.GetCommunicator();
    ModelPart::NodesContainerType& r_local_nodes = r_comm.LocalMesh().Nodes();

    int num_local_nodes = r_local_nodes.size();
    const int id_offset = 3*(r_comm.GetDataCommunicator().ScanSum(num_local_nodes) - num_local_nodes);

    int i = 1; // Dof Ids starting from one (to make sure that no dof with id 0 exists while error checking later)
    for (auto i_node = r_local_nodes.begin(); i_node != r_local_nodes.end(); ++i_node)
    {
        auto& r_dofs = i_node->GetDofs();
        for (auto i_dof = r_dofs.begin(); i_dof != r_dofs.end(); ++i_dof)
        {
            (*i_dof)->SetEquationId(id_offset + i);
            ++i;
        }
    }

    r_comm.SynchronizeDofs();
    for (auto i_node = r_model_part.NodesBegin(); i_node != r_model_part.NodesEnd(); ++i_node)
    {
        auto& r_dofs = i_node->GetDofs();
        for (auto i_dof = r_dofs.begin(); i_dof != r_dofs.end(); ++i_dof)
        {
            KRATOS_CHECK_NOT_EQUAL((*i_dof)->EquationId(), 0);
        }
    }
}


KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(ParallelFillCommunicatorExecution, KratosMPICoreFastSuite)
{
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("TestModelPart");
    r_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    MPIDataCommunicator comm_world(MPI_COMM_WORLD);
    Internals::ModelPartForMPICommunicatorTests(r_model_part, comm_world);

    auto& r_mpi_comm = r_model_part.GetCommunicator();
    unsigned int number_of_colors = r_mpi_comm.GetNumberOfColors();
    auto neighbor_indices = r_mpi_comm.NeighbourIndices();

    int neighbor;
    double local_index = comm_world.Rank(); // PARTITION_INDEX is a double
    for (unsigned int i = 0; i < number_of_colors; i++)
    {
        if ((neighbor = neighbor_indices[i]) > -1)
        {
            std::size_t interface_size = r_mpi_comm.InterfaceMeshes()[i].Nodes().size();
            std::size_t local_size = r_mpi_comm.LocalMeshes()[i].Nodes().size();
            std::size_t ghost_size = r_mpi_comm.GhostMeshes()[i].Nodes().size();
            KRATOS_CHECK_GREATER(interface_size, 0);
            KRATOS_CHECK_EQUAL(interface_size, local_size+ghost_size);
            double neighbor_index = neighbor; // PARTITION_INDEX is a double
            for (auto& node : r_mpi_comm.LocalMeshes()[i].Nodes())
            {
                KRATOS_CHECK_EQUAL(node.FastGetSolutionStepValue(PARTITION_INDEX,0), local_index);
            }
            for (auto& node : r_mpi_comm.GhostMeshes()[i].Nodes())
            {
                KRATOS_CHECK_EQUAL(node.FastGetSolutionStepValue(PARTITION_INDEX,0), neighbor_index);
            }
        }
    }
}

}
}