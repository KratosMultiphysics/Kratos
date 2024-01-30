//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <limits>

// External includes

/* Trilinos includes */
#include "Epetra_FEVector.h"

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "custom_strategies/convergencecriterias/trilinos_displacement_criteria.h"
#include "mpi/includes/mpi_data_communicator.h"
#include "mpi/utilities/parallel_fill_communicator.h"
#include "mpi/utilities/model_part_communicator_utilities.h"
#include "spaces/ublas_space.h"
#include "trilinos_space.h"

namespace Kratos::Testing
{
using GeometryType = Geometry<Node>;
using TrilinosSparseSpaceType = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>;
using TrilinosLocalSpaceType = UblasSpace<double, Matrix, Vector>;

using DofsArrayType = ModelPart::DofsArrayType;

using TrilinosDisplacementCriteriaType = TrilinosDisplacementCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType>;

void GenerateTestTrilinosDisplacementCriteriaModelPart(
    ModelPart& rModelPart,
    const unsigned int NumberOfNodes = 10
    )
{
    // The data communicator
    const DataCommunicator& r_comm = Testing::GetDefaultDataCommunicator();

    // Check communicator
    KRATOS_ERROR_IF_NOT(r_comm.IsDistributed()) << "Only distributed DataCommunicators can be used!" << std::endl;

    // Set MPI communicator
    ModelPartCommunicatorUtilities::SetMPICommunicator(rModelPart, r_comm);

    // Model part settings
    rModelPart.SetBufferSize(1);
    rModelPart.AddNodalSolutionStepVariable(PRESSURE);
    rModelPart.AddNodalSolutionStepVariable(PARTITION_INDEX);

    // MPI data
    const int rank = r_comm.Rank();
    const int size = r_comm.Size();

    // Create the auxiliary set of nodes
    for (unsigned int i_node = NumberOfNodes * rank; i_node < NumberOfNodes * (rank + 1); ++i_node) {
        auto p_node = rModelPart.CreateNewNode(i_node + 1, 0.0, 0.0, 0.0); // Coordinates do not matter in this test
        p_node->FastGetSolutionStepValue(PARTITION_INDEX) = rank;
        rModelPart.AddNode(p_node);
    }
    // Connected nodes
    if (rank < (size - 1)) {
        auto p_node = rModelPart.CreateNewNode((NumberOfNodes * (rank + 1)) + 1, 0.0, 0.0, 0.0); // Coordinates do not matter in this test
        p_node->FastGetSolutionStepValue(PARTITION_INDEX) = rank + 1;
        rModelPart.AddNode(p_node);
    }

    // Compute communication plan and fill communicator meshes correctly
    ParallelFillCommunicator(rModelPart, r_comm).Execute();
}

/**
 * Checks the displacement criteria
 */
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosDisplacementCriteria, KratosTrilinosApplicationMPITestSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("TestModelPart");

    // Set the test model part
    const unsigned int n_nodes = 10;
    GenerateTestTrilinosDisplacementCriteriaModelPart(r_model_part, n_nodes); // Create the geometry

    // MPI data
    const auto& r_comm = r_model_part.GetCommunicator().GetDataCommunicator();
    const int rank = r_comm.Rank();
    //const int size = r_comm.Size();

    // Add the corresponding test DOFs
    r_model_part.GetNodalSolutionStepVariablesList().AddDof(&PRESSURE);
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(PRESSURE);
    }

    // Create the displacement criteria
    const double rel_tol = 1.0e-3;
    const double abs_tol = 1.0e-5;
    auto displacement_criteria = TrilinosDisplacementCriteriaType(rel_tol, abs_tol);

    // Set the auxiliary arrays
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.pGetDof(PRESSURE)->SetEquationId(r_node.Id() - 1);
    }
    DofsArrayType aux_dof_set;
    aux_dof_set.reserve(r_model_part.NumberOfNodes());
    for (auto& r_node : r_model_part.Nodes()) {
        aux_dof_set.push_back(r_node.pGetDof(PRESSURE));
    }
    aux_dof_set.Sort();

    // Generate Epetra communicator
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    Epetra_MpiComm epetra_comm(raw_mpi_comm);

    // Generate map - use the "temp" array here
    const int temp_size = (n_nodes < 1000) ? 1000 : n_nodes;
    std::vector<int> temp_primary(temp_size, 0);
    for (IndexType i = 0; i != n_nodes; i++) {
        temp_primary[i] = (n_nodes * rank) + i;
    }
    Epetra_Map map(-1, n_nodes, temp_primary.data(), 0, epetra_comm);

    TrilinosSparseSpaceType::MatrixPointerType pA;
    TrilinosSparseSpaceType::MatrixType& rA  = *pA;  /// The LHS matrix of the system of equations (Only required to match the API)
    TrilinosSparseSpaceType::VectorType Dx(map);     /// The increment in the solution
    TrilinosSparseSpaceType::VectorPointerType pb;
    TrilinosSparseSpaceType::VectorType& rb  = *pb;  /// The RHS vector of the system of equations (Only required to match the API)

    // Set the auxiliary fake data to check the convergence for
    const double aux_constant = 1.0e-3;

    // Set the auxiliary fake data to check the convergence for (initialization)
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.FastGetSolutionStepValue(PRESSURE) = aux_constant;
        if (r_node.FastGetSolutionStepValue(PARTITION_INDEX) == rank) {
            const IndexType index = r_node.Id() - 1;
            const double aux_val = aux_constant;
            TrilinosSparseSpaceType::SetValue(Dx, index, aux_val);
        }
    }

    // Initialize the solution step
    displacement_criteria.InitializeSolutionStep(r_model_part, aux_dof_set, rA, Dx, rb);

    // Check convergence (failing)
    bool convergence = displacement_criteria.PostCriteria(r_model_part, aux_dof_set, rA, Dx, rb);
    KRATOS_EXPECT_FALSE(convergence)

    // Set the auxiliary fake data to check the convergence for (passing)
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.FastGetSolutionStepValue(PRESSURE) = aux_constant;
        if (r_node.FastGetSolutionStepValue(PARTITION_INDEX) == rank) {
            const IndexType index = r_node.Id() - 1;
            const double aux_val = aux_constant/1000.0;
            TrilinosSparseSpaceType::SetValue(Dx, index, aux_val);
        }
    }

    // Check convergence (passing)
    convergence = displacement_criteria.PostCriteria(r_model_part, aux_dof_set, rA, Dx, rb);
    KRATOS_EXPECT_TRUE(convergence)
}

} // namespace Kratos::Testing
