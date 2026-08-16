//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || "__"| | | | "_ \\ / _ \\/ __|
//           | || |  | | | | | | | (_) \\__
//           |_||_|  |_|_|_|_| |_|\\___/|___/ APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <limits>

// External includes

// Project includes
#include "tests/cpp_tests/trilinos_fast_suite.h"
#include "containers/model.h"
#include "custom_strategies/convergencecriterias/trilinos_residual_criteria.h"
#include "mpi/includes/mpi_data_communicator.h"
#include "mpi/utilities/parallel_fill_communicator.h"
#include "mpi/utilities/model_part_communicator_utilities.h"
#include "spaces/ublas_space.h"
#include "trilinos_space_experimental.h"

namespace Kratos::Testing
{
using GeometryType = Geometry<Node>;
using TrilinosSparseSpaceType = TrilinosSpaceExperimental<Tpetra::FECrsMatrix<>, Tpetra::FEMultiVector<>>;
using TrilinosLocalSpaceType = UblasSpace<double, Matrix, Vector>;

using DofsArrayType = ModelPart::DofsArrayType;

using TrilinosResidualCriteriaType = TrilinosResidualCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType>;

void GenerateTestTrilinosExperimentalResidualCriteriaModelPart(
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
 * Checks the residual criteria
 */
KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalResidualCriteria, KratosTrilinosApplicationMPITestSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("TestModelPart");

    // Set the test model part
    const unsigned int n_nodes = 10;
    GenerateTestTrilinosExperimentalResidualCriteriaModelPart(r_model_part, n_nodes); // Create the geometry

    // MPI data
    const auto& r_comm = r_model_part.GetCommunicator().GetDataCommunicator();
    const int rank = r_comm.Rank();

    // Add the corresponding test DOFs
    r_model_part.GetNodalSolutionStepVariablesList().AddDof(&PRESSURE);
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(PRESSURE);
    }

    // Create the residual criteria
    const double rel_tol = 1.0e-3;
    const double abs_tol = 1.0e-5;
    auto residual_criteria = TrilinosResidualCriteriaType(rel_tol, abs_tol);

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

    // Tpetra communicator
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    Teuchos::RCP<const Teuchos::MpiComm<int>> p_tpetra_comm = Teuchos::rcp(new Teuchos::MpiComm<int>(raw_mpi_comm));

    // Generate map
    using LO = TrilinosSparseSpaceType::LO;
    using GO = TrilinosSparseSpaceType::GO;
    using NT = TrilinosSparseSpaceType::NT;
    std::vector<GO> indices(n_nodes);
    for (IndexType i = 0; i != n_nodes; i++) {
        indices[i] = (n_nodes * rank) + i;
    }
    Teuchos::RCP<const Tpetra::Map<LO, GO, NT>> p_map = Teuchos::rcp(new Tpetra::Map<LO, GO, NT>(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), Teuchos::ArrayView<const GO>(indices), 0, p_tpetra_comm));

    // Not used in this criteria but needed for API
    auto p_map_aux = Teuchos::rcp(new TrilinosSparseSpaceType::MapType(0, 0, p_tpetra_comm));
    auto p_graph_aux = Teuchos::rcp(new TrilinosSparseSpaceType::GraphType(p_map_aux, p_map_aux, 0));
    p_graph_aux->fillComplete();
    TrilinosSparseSpaceType::MatrixType rA(p_graph_aux);
    auto p_Dx = TrilinosSparseSpaceType::CreateVector(p_map_aux); // Not used in this criteria but needed for API
    auto& Dx = *p_Dx;
    auto p_rb = TrilinosSparseSpaceType::CreateVector(p_map);
    auto& rb = *p_rb;

    // Set the auxiliary fake data to check the convergence for
    const double aux_constant = 1.0e-3;

    // Set the auxiliary fake data to check the convergence for (initialization)
    for (auto& r_node : r_model_part.Nodes()) {
        if (r_node.FastGetSolutionStepValue(PARTITION_INDEX) == rank) {
            const IndexType index = r_node.Id() - 1;
            const double aux_val = aux_constant;
            TrilinosSparseSpaceType::SetValue(rb, index, aux_val);
        }
    }

    // Initialize the solution step
    residual_criteria.InitializeSolutionStep(r_model_part, aux_dof_set, rA, Dx, rb);

    // Check convergence (failing)
    bool convergence = residual_criteria.PostCriteria(r_model_part, aux_dof_set, rA, Dx, rb);
    KRATOS_EXPECT_FALSE(convergence)

    // Set the auxiliary fake data to check the convergence for (passing)
    for (auto& r_node : r_model_part.Nodes()) {
        if (r_node.FastGetSolutionStepValue(PARTITION_INDEX) == rank) {
            const IndexType index = r_node.Id() - 1;
            const double aux_val = aux_constant/1000.0;
            TrilinosSparseSpaceType::SetValue(rb, index, aux_val);
        }
    }

    // Check convergence (passing)
    convergence = residual_criteria.PostCriteria(r_model_part, aux_dof_set, rA, Dx, rb);
    KRATOS_EXPECT_TRUE(convergence)
}

} // namespace Kratos::Testing
