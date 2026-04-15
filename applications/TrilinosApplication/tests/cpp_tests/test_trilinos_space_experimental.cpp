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

// External includes

// Project includes
#include "trilinos_space_experimental.h"
#include "tests/cpp_tests/trilinos_fast_suite.h"
#include "containers/model.h"
#include "mpi/includes/mpi_data_communicator.h"
#include "trilinos_cpp_test_experimental_utilities.h"
#include "utilities/math_utils.h"

namespace Kratos::Testing
{

/// Basic definitions
using TrilinosSparseMatrix = Tpetra::FECrsMatrix<>;
using TrilinosSparseVector = Tpetra::FEMultiVector<>;
using TrilinosSparseSpaceType = TrilinosSpaceExperimental<TrilinosSparseMatrix, TrilinosSparseVector>;
using TrilinosLocalSpaceType = UblasSpace<double, Matrix, Vector>;

using TrilinosSparseMatrixType = TrilinosSparseSpaceType::MatrixType;
using TrilinosVectorType = TrilinosSparseSpaceType::VectorType;

using TrilinosMatrixPointerType = TrilinosSparseSpaceType::MatrixPointerType;
using TrilinosVectorPointerType = TrilinosSparseSpaceType::VectorPointerType;

using TrilinosLocalMatrixType = TrilinosLocalSpaceType::MatrixType;
using TrilinosLocalVectorType = TrilinosLocalSpaceType::VectorType;

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalSizeVector, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto vector = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size);

    // Check
    KRATOS_EXPECT_EQ(static_cast<std::size_t>(size), TrilinosSparseSpaceType::Size(*vector));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalSizeMatrix, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size);

    // Check
    KRATOS_EXPECT_EQ(static_cast<std::size_t>(size), TrilinosSparseSpaceType::Size1(*matrix));
    KRATOS_EXPECT_EQ(static_cast<std::size_t>(size), TrilinosSparseSpaceType::Size2(*matrix));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalDotProduct, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto vector1 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size);
    auto vector2 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size);
    double ref = 0.0;
    for (int i = 0; i < size; ++i) {
        ref += std::pow(static_cast<double>(i), 2);
    }

    // Check
    KRATOS_EXPECT_DOUBLE_EQ(ref, TrilinosSparseSpaceType::Dot(*vector1, *vector2));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalMaxMin, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto vector = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size);

    // Check
    KRATOS_EXPECT_DOUBLE_EQ(0.0, TrilinosSparseSpaceType::Min(*vector));
    KRATOS_EXPECT_DOUBLE_EQ(static_cast<double>(size - 1), TrilinosSparseSpaceType::Max(*vector));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalTwoNormVector, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto vector = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size);
    double ref = 0.0;
    for (int i = 0; i < size; ++i) {
        ref += std::pow(static_cast<double>(i), 2);
    }
    ref = std::sqrt(ref);

    // Check
    KRATOS_EXPECT_DOUBLE_EQ(ref, TrilinosSparseSpaceType::TwoNorm(*vector));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalTwoNormMatrix1, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 2 * r_comm.Size();
    auto matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size);
    double ref = 0.0;
    for (int i = 0; i < size; ++i) {
        ref += std::pow(static_cast<double>(i), 2);
    }
    ref = std::sqrt(ref);

    // Check
    KRATOS_EXPECT_DOUBLE_EQ(ref, TrilinosSparseSpaceType::TwoNorm(*matrix));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalTwoNormMatrix2, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 2 * r_comm.Size();
    auto matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 10.0, true, 100);
    auto local_matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalMatrix(size, 10.0, true);

    // Check
    KRATOS_EXPECT_DOUBLE_EQ(TrilinosLocalSpaceType::TwoNorm(local_matrix), TrilinosSparseSpaceType::TwoNorm(*matrix));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalMultMatrixVector, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 2 * r_comm.Size();
    auto matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 10.0, true, 100);
    auto local_matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalMatrix(size, 10.0, true);
    auto vector = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size, 0.0);
    auto local_vector = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalVector(size, 0.0);

    // Tpetra communicator
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    TrilinosSparseSpaceType::CommunicatorPointerType tpetra_comm = Teuchos::rcp(new TrilinosSparseSpaceType::CommunicatorType(raw_mpi_comm));

    // Create a map
    const int global_elems = size;
    TrilinosSparseSpaceType::MapPointerType map = Teuchos::rcp(new TrilinosSparseSpaceType::MapType(global_elems, 0, tpetra_comm));

    // Create an Tpetra_Vector
    auto p_mult = TrilinosSparseSpaceType::CreateVector(map);
    auto& mult = *p_mult;

    // Solution
    TrilinosSparseSpaceType::Mult(*matrix, *vector, mult);

    // Check
    const TrilinosLocalVectorType multiply_reference = prod(local_matrix, local_vector);
    TrilinosCPPTestExperimentalUtilities::CheckSparseVectorFromLocalVector(mult, multiply_reference);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalTransposeMultMatrixMatrix, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 2 * r_comm.Size();
    auto matrix_1 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 10.0, true, 100);
    auto local_matrix_1 = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalMatrix(size, 10.0, true);
    auto matrix_2 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 20.0, true, 100);
    auto local_matrix_2 = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalMatrix(size, 20.0, true);

        // Create an Tpetra_Matrix
    // Create an Tpetra_Matrix with enough capacity for product
    using LO = TrilinosSparseSpaceType::LO;
    using GO = TrilinosSparseSpaceType::GO;
    Teuchos::RCP<TrilinosSparseSpaceType::GraphType> p_graph_mult = Teuchos::rcp(new TrilinosSparseSpaceType::GraphType(matrix_1->getRowMap(), matrix_1->getRowMap(), size));
    for (LO i = 0; i < static_cast<LO>(matrix_1->getNodeNumRows()); ++i) {
        const GO global_row = matrix_1->getRowMap()->getGlobalElement(i);
        std::vector<GO> indices(size);
        for (int j = 0; j < size; ++j) indices[j] = j;
        p_graph_mult->insertGlobalIndices(global_row, Teuchos::ArrayView<const GO>(indices));
    }
    p_graph_mult->fillComplete();
    TrilinosSparseMatrixType mult(p_graph_mult);

    // Solution
    TrilinosSparseSpaceType::TransposeMult(*matrix_1, *matrix_2, mult, {true, false});

    // Check
    const TrilinosLocalMatrixType multiply_reference = prod(trans(local_matrix_1), local_matrix_2);
    TrilinosCPPTestExperimentalUtilities::CheckSparseMatrixFromLocalMatrix(mult, multiply_reference);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalBtDBProductOperation, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 2 * r_comm.Size();
    auto matrix_1 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 10.0, true, 100);
    auto local_matrix_1 = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalMatrix(size, 10.0, true);
    auto matrix_2 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 20.0, true, 100);
    auto local_matrix_2 = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalMatrix(size, 20.0, true);

    // Create an Tpetra_Matrix with correct sparsity for product
    // Create an Tpetra_Matrix with enough capacity for product
    using LO = TrilinosSparseSpaceType::LO;
    using GO = TrilinosSparseSpaceType::GO;
    Teuchos::RCP<TrilinosSparseSpaceType::GraphType> p_graph_mult = Teuchos::rcp(new TrilinosSparseSpaceType::GraphType(matrix_1->getRowMap(), matrix_1->getRowMap(), size));
    for (LO i = 0; i < static_cast<LO>(matrix_1->getNodeNumRows()); ++i) {
        const GO global_row = matrix_1->getRowMap()->getGlobalElement(i);
        std::vector<GO> indices(size);
        for (int j = 0; j < size; ++j) indices[j] = j;
        p_graph_mult->insertGlobalIndices(global_row, Teuchos::ArrayView<const GO>(indices));
    }
    p_graph_mult->fillComplete();
    TrilinosSparseMatrixType mult(p_graph_mult);

    // Solution
    TrilinosSparseSpaceType::BtDBProductOperation(mult, *matrix_1, *matrix_2);

    // Check
    TrilinosLocalMatrixType multiply_reference;
    MathUtils<double>::BtDBProductOperation(multiply_reference, local_matrix_1, local_matrix_2);
    TrilinosCPPTestExperimentalUtilities::CheckSparseMatrixFromLocalMatrix(mult, multiply_reference);

    // Non zero matrix
    Teuchos::RCP<TrilinosSparseSpaceType::GraphType> graph_second = Teuchos::rcp(new TrilinosSparseSpaceType::GraphType(matrix_1->getRowMap(), matrix_1->getRowMap(), size));
    for (LO i = 0; i < static_cast<LO>(matrix_1->getNodeNumRows()); ++i) {
        const GO global_row = matrix_1->getRowMap()->getGlobalElement(i);
        std::vector<GO> indices(size);
        for (int j = 0; j < size; ++j) indices[j] = j;
        graph_second->insertGlobalIndices(global_row, Teuchos::ArrayView<const GO>(indices));
    }
    graph_second->fillComplete();
    auto second_mult = Teuchos::rcp(new TrilinosSparseMatrixType(graph_second));

    // Solution
    TrilinosSparseSpaceType::BtDBProductOperation(*second_mult, *matrix_1, *matrix_2);
    TrilinosCPPTestExperimentalUtilities::CheckSparseMatrixFromLocalMatrix(*second_mult, multiply_reference);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalBDBtProductOperation, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 2 * r_comm.Size();
    auto matrix_1 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 10.0, true, 100);
    auto local_matrix_1 = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalMatrix(size, 10.0, true);
    auto matrix_2 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 20.0, true, 100);
    auto local_matrix_2 = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalMatrix(size, 20.0, true);

    // Create an Tpetra_Matrix with correct sparsity for product
    // Create an Tpetra_Matrix with enough capacity for product
    using LO = TrilinosSparseSpaceType::LO;
    using GO = TrilinosSparseSpaceType::GO;
    Teuchos::RCP<TrilinosSparseSpaceType::GraphType> p_graph_mult = Teuchos::rcp(new TrilinosSparseSpaceType::GraphType(matrix_1->getRowMap(), matrix_1->getRowMap(), size));
    for (LO i = 0; i < static_cast<LO>(matrix_1->getNodeNumRows()); ++i) {
        const GO global_row = matrix_1->getRowMap()->getGlobalElement(i);
        std::vector<GO> indices(size);
        for (int j = 0; j < size; ++j) indices[j] = j;
        p_graph_mult->insertGlobalIndices(global_row, Teuchos::ArrayView<const GO>(indices));
    }
    p_graph_mult->fillComplete();
    TrilinosSparseMatrixType mult(p_graph_mult);

    // Solution
    TrilinosSparseSpaceType::BDBtProductOperation(mult, *matrix_1, *matrix_2);

    // Check
    TrilinosLocalMatrixType multiply_reference;
    MathUtils<double>::BDBtProductOperation(multiply_reference, local_matrix_1, local_matrix_2);
    TrilinosCPPTestExperimentalUtilities::CheckSparseMatrixFromLocalMatrix(mult, multiply_reference);

    // Non zero matrix
    Teuchos::RCP<TrilinosSparseSpaceType::GraphType> graph_second = Teuchos::rcp(new TrilinosSparseSpaceType::GraphType(matrix_1->getRowMap(), matrix_1->getRowMap(), size));
    for (LO i = 0; i < static_cast<LO>(matrix_1->getNodeNumRows()); ++i) {
        const GO global_row = matrix_1->getRowMap()->getGlobalElement(i);
        std::vector<GO> indices(size);
        for (int j = 0; j < size; ++j) indices[j] = j;
        graph_second->insertGlobalIndices(global_row, Teuchos::ArrayView<const GO>(indices));
    }
    graph_second->fillComplete();
    auto second_mult = Teuchos::rcp(new TrilinosSparseMatrixType(graph_second));

    // Solution
    TrilinosSparseSpaceType::BDBtProductOperation(*second_mult, *matrix_1, *matrix_2);
    TrilinosCPPTestExperimentalUtilities::CheckSparseMatrixFromLocalMatrix(*second_mult, multiply_reference);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalInplaceMult, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto vector = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size);
    auto local_vector = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalVector(size);

    // Multiply
    const double mult = 2.0;

    // Solution
    TrilinosSparseSpaceType::InplaceMult(*vector, mult);

    // Check
    local_vector *= mult;
    TrilinosCPPTestExperimentalUtilities::CheckSparseVectorFromLocalVector(*vector, local_vector);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalAssign, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto vector_1 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size);
    auto vector_2 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size, 1.0);
    auto local_vector_1 = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalVector(size);
    auto local_vector_2 = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalVector(size, 1.0);

    // Multiply
    const double mult = 2.0;

    // Solution
    TrilinosSparseSpaceType::Assign(*vector_1, mult, *vector_2);

    // Check
    local_vector_1 = mult * local_vector_2;
    TrilinosCPPTestExperimentalUtilities::CheckSparseVectorFromLocalVector(*vector_1, local_vector_1);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalUnaliasedAdd, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto vector_1 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size);
    auto vector_2 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size, 1.0);
    auto local_vector_1 = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalVector(size);
    auto local_vector_2 = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalVector(size, 1.0);

    // Multiply
    const double mult = 2.0;

    // Solution
    TrilinosSparseSpaceType::UnaliasedAdd(*vector_1, mult, *vector_2);

    // Check
    local_vector_1 += mult * local_vector_2;
    TrilinosCPPTestExperimentalUtilities::CheckSparseVectorFromLocalVector(*vector_1, local_vector_1);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalScaleAndAdd1, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto vector_1 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size);
    auto vector_2 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size, 1.0);
    auto vector_3 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size, 2.0);
    auto local_vector_1 = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalVector(size);
    auto local_vector_2 = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalVector(size, 1.0);
    auto local_vector_3 = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalVector(size, 2.0);

    // Multiply
    const double mult_1 = 2.0;
    const double mult_2 = 1.5;

    // Solution
    TrilinosSparseSpaceType::ScaleAndAdd(mult_1, *vector_2, mult_2, *vector_3, *vector_1);

    // Check
    local_vector_1 = mult_1 * local_vector_2 + mult_2 * local_vector_3;
    TrilinosCPPTestExperimentalUtilities::CheckSparseVectorFromLocalVector(*vector_1, local_vector_1);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalScaleAndAdd2, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto vector_1 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size);
    auto vector_2 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size, 1.0);
    auto local_vector_1 = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalVector(size);
    auto local_vector_2 = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalVector(size, 1.0);

    // Multiply
    const double mult_1 = 2.0;
    const double mult_2 = 1.5;

    // Solution
    TrilinosSparseSpaceType::ScaleAndAdd(mult_1, *vector_2, mult_2, *vector_1);

    // Check
    local_vector_1 = mult_1 * local_vector_2 + mult_2 * local_vector_1;
    TrilinosCPPTestExperimentalUtilities::CheckSparseVectorFromLocalVector(*vector_1, local_vector_1);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalSet, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto vector = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size);
    auto local_vector = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalVector(size);

    // Set
    const double value = 2.0;

    // Solution
    TrilinosSparseSpaceType::Set(*vector, value);

    // Check
    for (int i = 0; i < size; ++i) local_vector[i] = value;
    TrilinosCPPTestExperimentalUtilities::CheckSparseVectorFromLocalVector(*vector, local_vector);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalSetToZeroMatrix, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size);
    TrilinosLocalMatrixType local_matrix = ZeroMatrix(size, size);

    // Solution
    TrilinosSparseSpaceType::SetToZero(*matrix);

    // Check
    TrilinosCPPTestExperimentalUtilities::CheckSparseMatrixFromLocalMatrix(*matrix, local_matrix);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalSetToZeroVector, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto vector = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size);
    TrilinosLocalVectorType local_vector = ZeroVector(size);

    // Solution
    TrilinosSparseSpaceType::SetToZero(*vector);

    // Check
    TrilinosCPPTestExperimentalUtilities::CheckSparseVectorFromLocalVector(*vector, local_vector);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalCopyMatrixValues, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto matrix_1 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 10.0, true, 100);
    auto matrix_2 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 0.0);
    auto local_matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalMatrix(size, 0.0);

    // Solution
    TrilinosSparseSpaceType::CopyMatrixValues(*matrix_1, *matrix_2);

    // Check
    TrilinosCPPTestExperimentalUtilities::CheckSparseMatrixFromLocalMatrix(*matrix_1, local_matrix);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalCombineMatricesGraphs, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto matrix_1 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 0.0);
    auto matrix_2 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 10.0, true, 100);
    auto local_matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalMatrix(size, 10.0, true);

    // Creating new matrix
    const auto combined_graph = TrilinosSparseSpaceType::CombineMatricesGraphs(*matrix_1, *matrix_2);
    TrilinosMatrixPointerType copied_matrix = Teuchos::rcp(new TrilinosSparseMatrixType(combined_graph));

    // Solution
    TrilinosSparseSpaceType::CopyMatrixValues(*copied_matrix, *matrix_2);

    // Check
    TrilinosCPPTestExperimentalUtilities::CheckSparseMatrixFromLocalMatrix(*copied_matrix, local_matrix);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalCheckAndCorrectZeroDiagonalValues, KratosTrilinosApplicationMPITestSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");
    auto& r_process_info = r_model_part.GetProcessInfo();
    r_process_info.SetValue(BUILD_SCALE_FACTOR, 1.0);

    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 12;
    auto matrix12x12 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size);

    // Generate Tpetra communicator
    KRATOS_ERROR_IF_NOT(r_comm.IsDistributed()) << "Only distributed DataCommunicators can be used!" << std::endl;
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    TrilinosSparseSpaceType::CommunicatorPointerType tpetra_comm = Teuchos::rcp(new TrilinosSparseSpaceType::CommunicatorType(raw_mpi_comm));

    // Dummy vector
    auto map = Teuchos::rcp(new TrilinosSparseSpaceType::MapType(size, size, tpetra_comm));
    auto p_vector12 = TrilinosSparseSpaceType::CreateVector(map);
    auto& vector12 = *p_vector12;

    // Test the norm of the matrix
    double norm = 0.0;
    norm = TrilinosSparseSpaceType::CheckAndCorrectZeroDiagonalValues(r_process_info, *matrix12x12, vector12, SCALING_DIAGONAL::NO_SCALING);
    KRATOS_EXPECT_DOUBLE_EQ(norm, 1.0);
    if (r_comm.Rank() == 0) {
        auto localMatrix = matrix12x12->getLocalMatrixHost();
        auto valuesView = localMatrix.values;
        KRATOS_EXPECT_DOUBLE_EQ(valuesView(0), 1.0);
    }
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalIsDistributed, KratosTrilinosApplicationMPITestSuite)
{
    KRATOS_EXPECT_TRUE(TrilinosSparseSpaceType::IsDistributed());
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalGetScaleNorm, KratosTrilinosApplicationMPITestSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");
    auto& r_process_info = r_model_part.GetProcessInfo();
    r_process_info.SetValue(BUILD_SCALE_FACTOR, 3.0);

    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 12;
    auto matrix12x12 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 1.0);

    // Test the norm of the matrix
    double norm = 0.0;
    norm = TrilinosSparseSpaceType::GetScaleNorm(r_process_info, *matrix12x12, SCALING_DIAGONAL::NO_SCALING);
    KRATOS_EXPECT_DOUBLE_EQ(norm, 1.0);
    norm = TrilinosSparseSpaceType::GetScaleNorm(r_process_info, *matrix12x12, SCALING_DIAGONAL::CONSIDER_PRESCRIBED_DIAGONAL);
    KRATOS_EXPECT_DOUBLE_EQ(norm, 3.0);
    norm = TrilinosSparseSpaceType::GetScaleNorm(r_process_info, *matrix12x12, SCALING_DIAGONAL::CONSIDER_NORM_DIAGONAL);
    KRATOS_EXPECT_NEAR(norm, 2.124591464, 1.0e-6);
    norm = TrilinosSparseSpaceType::GetScaleNorm(r_process_info, *matrix12x12, SCALING_DIAGONAL::CONSIDER_MAX_DIAGONAL);
    KRATOS_EXPECT_DOUBLE_EQ(norm, 12.0);
    norm = TrilinosSparseSpaceType::GetAveragevalueDiagonal(*matrix12x12);
    KRATOS_EXPECT_DOUBLE_EQ(norm, 6.5);
    norm = TrilinosSparseSpaceType::GetMinDiagonal(*matrix12x12);
    KRATOS_EXPECT_DOUBLE_EQ(norm, 1.0);
}


KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalScaleAndAddMatrix, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 12;
    auto matrix_1 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 10.0, true, 100);
    auto matrix_2 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 20.0, true, 100);

    if (r_comm.Rank() == 0) {
        auto v1 = matrix_1->getLocalMatrixHost();
        auto v2 = matrix_2->getLocalMatrixHost();
        KRATOS_WATCH(v1.values(0));
        KRATOS_WATCH(v2.values(0));
    }

    // Solution
    TrilinosSparseSpaceType::ScaleAndAdd(2.0, *matrix_1, 3.0, *matrix_2);

    // Check
    auto local_matrix_1 = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalMatrix(size, 10.0, true);
    auto local_matrix_2 = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalMatrix(size, 20.0, true);
    TrilinosLocalMatrixType local_reference = 2.0 * local_matrix_1 + 3.0 * local_matrix_2;

    TrilinosCPPTestExperimentalUtilities::CheckSparseMatrixFromLocalMatrix(*matrix_2, local_reference);
}

/* KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalMatrixMarket, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 12;
    auto matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 10.0, true, 100);

    // Write to file
    const std::string file_name = "test_matrix.mm";
    TrilinosSparseSpaceType::WriteMatrixMarketMatrix(file_name.c_str(), *matrix, false);

    // Read from file
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    Teuchos::MpiComm<int> tpetra_comm(raw_mpi_comm);
    auto read_matrix = TrilinosSparseSpaceType::ReadMatrixMarket(file_name, tpetra_comm);

    // Check
    auto local_matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalMatrix(size, 10.0, true);
    TrilinosCPPTestExperimentalUtilities::CheckSparseMatrixFromLocalMatrix(*read_matrix, local_matrix);

    // Clean up
    if (r_comm.Rank() == 0) {
        std::remove(file_name.c_str());
    }
} */


/* KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalMatrixMarketVector, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 12;
    auto vector = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size);

    // Write to file
    const std::string file_name = "test_vector.mm";
    TrilinosSparseSpaceType::WriteMatrixMarketVector(file_name.c_str(), *vector);

    // Tpetra communicator
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    TrilinosSparseSpaceType::CommunicatorPointerType tpetra_comm = Teuchos::rcp(new TrilinosSparseSpaceType::CommunicatorType(raw_mpi_comm));

    // Read from file
    auto read_vector = TrilinosSparseSpaceType::ReadMatrixMarketVector(file_name, tpetra_comm, size);

    // Check
    auto local_vector = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalVector(size);
    TrilinosCPPTestExperimentalUtilities::CheckSparseVectorFromLocalVector(*read_vector, local_vector);

    // Clean up
    if (r_comm.Rank() == 0) {
        std::remove(file_name.c_str());
    }
} */


KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalTransposeMultMatrixVector, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 2 * r_comm.Size();
    auto matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 10.0, true, 100);
    auto local_matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalMatrix(size, 10.0, true);
    auto vector = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size, 10.0);
    auto local_vector = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalVector(size, 10.0);

    // Create an Tpetra_Vector
    auto p_mult = TrilinosSparseSpaceType::CreateVector(matrix->getDomainMap());
    auto& mult = *p_mult;

    // Solution
    TrilinosSparseSpaceType::TransposeMult(*matrix, *vector, mult);

    // Check
    const TrilinosLocalVectorType multiply_reference = prod(trans(local_matrix), local_vector);
    TrilinosCPPTestExperimentalUtilities::CheckSparseVectorFromLocalVector(mult, multiply_reference);
}


KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalMultMatrixMatrix, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 2 * r_comm.Size();
    auto matrix_1 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 10.0, true, 100);
    auto local_matrix_1 = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalMatrix(size, 10.0, true);
    auto matrix_2 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 20.0, true, 100);
    auto local_matrix_2 = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalMatrix(size, 20.0, true);

    // Create an Tpetra_Matrix with enough capacity for product
    using LO = TrilinosSparseSpaceType::LO;
    using GO = TrilinosSparseSpaceType::GO;
    Teuchos::RCP<TrilinosSparseSpaceType::GraphType> p_graph_mult = Teuchos::rcp(new TrilinosSparseSpaceType::GraphType(matrix_1->getRowMap(), matrix_1->getRowMap(), size));
    for (LO i = 0; i < static_cast<LO>(matrix_1->getNodeNumRows()); ++i) {
        const GO global_row = matrix_1->getRowMap()->getGlobalElement(i);
        std::vector<GO> indices(size);
        for (int j = 0; j < size; ++j) indices[j] = j;
        p_graph_mult->insertGlobalIndices(global_row, Teuchos::ArrayView<const GO>(indices));
    }
    p_graph_mult->fillComplete();
    TrilinosSparseMatrixType mult(p_graph_mult);

    // Solution
    TrilinosSparseSpaceType::Mult(*matrix_1, *matrix_2, mult);

    // Check
    const TrilinosLocalMatrixType multiply_reference = prod(local_matrix_1, local_matrix_2);
    TrilinosCPPTestExperimentalUtilities::CheckSparseMatrixFromLocalMatrix(mult, multiply_reference);
}


KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalIsDistributedSpace, KratosTrilinosApplicationMPITestSuite)
{
    KRATOS_EXPECT_TRUE(TrilinosSparseSpaceType::IsDistributedSpace());
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalLinearAlgebraLibrary, KratosTrilinosApplicationMPITestSuite)
{
    KRATOS_EXPECT_EQ(TrilinosSparseSpaceType::LinearAlgebraLibrary(), TrilinosLinearAlgebraLibrary::TPETRA);
}



KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalBtDBProductOperationRealCase, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 6;

    // Generate A matrix
    std::vector<int> row_indexes = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5};
    std::vector<int> column_indexes = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 2, 3, 4, 5, 2, 3, 4, 5};
    std::vector<double> values = {2069000000.0, 0.0, -2069000000.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2069000000.0, 0.0, 4138000000.0, 0.0, -2069000000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2069000000.0, 0.0, 2069000000.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    auto A = TrilinosCPPTestExperimentalUtilities::GenerateSparseMatrix(r_comm, size, row_indexes, column_indexes, values);

    // Generate T matrix
    row_indexes = {0,1,2,3,4,4,5};
    column_indexes = {0,1,2,3,2,4,5};
    values = {1,1,1,1,1,0,1};
    auto T = TrilinosCPPTestExperimentalUtilities::GenerateSparseMatrix(r_comm, size, row_indexes, column_indexes, values);

    /* Intermediate multiplication */

    // Create an Tpetra_Matrix with enough capacity
    using LO = TrilinosSparseSpaceType::LO;
    using GO = TrilinosSparseSpaceType::GO;
    Teuchos::RCP<TrilinosSparseSpaceType::GraphType> p_graph_aux = Teuchos::rcp(new TrilinosSparseSpaceType::GraphType(A->getRowMap(), A->getRowMap(), size));
    for (LO i = 0; i < static_cast<LO>(A->getNodeNumRows()); ++i) {
        const GO global_row = A->getRowMap()->getGlobalElement(i);
        std::vector<GO> indices(size);
        for (int j = 0; j < size; ++j) indices[j] = j;
        p_graph_aux->insertGlobalIndices(global_row, Teuchos::ArrayView<const GO>(indices));
    }
    p_graph_aux->fillComplete();
    TrilinosSparseMatrixType aux(p_graph_aux);

    // First multiplication (aux = T^T * A)
    TrilinosSparseSpaceType::TransposeMult(*T, *A, aux, {true, false});

    // Values to check
    row_indexes = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5};
    column_indexes = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 2, 3, 4, 5, 2, 3, 4, 5};
    values = {2069000000.0, 0.0, -2069000000.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2069000000.0, 0.0, 2069000000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // Check
    TrilinosCPPTestExperimentalUtilities::CheckSparseMatrix(aux, row_indexes, column_indexes, values);

    // Compute T^T A T
    // We create a result matrix with enough capacity
    Teuchos::RCP<TrilinosSparseSpaceType::GraphType> p_graph_res = Teuchos::rcp(new TrilinosSparseSpaceType::GraphType(A->getRowMap(), A->getRowMap(), size));
    for (LO i = 0; i < static_cast<LO>(A->getNodeNumRows()); ++i) {
        const GO global_row = A->getRowMap()->getGlobalElement(i);
        std::vector<GO> indices(size);
        for (int j = 0; j < size; ++j) indices[j] = j;
        p_graph_res->insertGlobalIndices(global_row, Teuchos::ArrayView<const GO>(indices));
    }
    p_graph_res->fillComplete();
    TrilinosSparseMatrixType res(p_graph_res);

    TrilinosSparseSpaceType::BtDBProductOperation(res, *A, *T);

    // Values to check
    row_indexes = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5};
    column_indexes = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 2, 3, 4, 5, 2, 3, 4, 5};
    values = {2069000000.0, 0.0, -2069000000.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2069000000.0, 0.0, 2069000000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // Check
    TrilinosCPPTestExperimentalUtilities::CheckSparseMatrix(res, row_indexes, column_indexes, values);
}


KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalGetRank, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    Teuchos::MpiComm<int> tpetra_comm(raw_mpi_comm);
    KRATOS_EXPECT_EQ(r_comm.Rank(), TrilinosSparseSpaceType::GetRank(tpetra_comm));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalIsNull, KratosTrilinosApplicationMPITestSuite)
{
    TrilinosSparseSpaceType::MatrixPointerType pA = Teuchos::null;
    KRATOS_EXPECT_TRUE(TrilinosSparseSpaceType::IsNull(pA));
    pA = TrilinosSparseSpaceType::CreateEmptyMatrixPointer();
    KRATOS_EXPECT_TRUE(TrilinosSparseSpaceType::IsNull(pA));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalCreateEmptyPointer, KratosTrilinosApplicationMPITestSuite)
{
    auto p_map = TrilinosSparseSpaceType::CreateEmptyMapPointer();
    KRATOS_EXPECT_TRUE(TrilinosSparseSpaceType::IsNull(p_map));
    auto p_A = TrilinosSparseSpaceType::CreateEmptyMatrixPointer();
    KRATOS_EXPECT_TRUE(TrilinosSparseSpaceType::IsNull(p_A));
    auto p_v = TrilinosSparseSpaceType::CreateEmptyVectorPointer();
    KRATOS_EXPECT_TRUE(TrilinosSparseSpaceType::IsNull(p_v));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalGetMap, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    auto vector = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size);
    auto& r_map = TrilinosSparseSpaceType::GetMap(*vector);
    KRATOS_EXPECT_EQ(static_cast<std::size_t>(size), r_map.getGlobalNumElements());
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalGetCommunicator, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    auto vector = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size);
    auto& r_trilinos_comm = TrilinosSparseSpaceType::GetCommunicator(*vector);
    KRATOS_EXPECT_EQ(r_comm.Rank(), r_trilinos_comm.getRank());
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalGlobalAssemble, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    auto vector = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size);
    auto matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size);
    TrilinosSparseSpaceType::GlobalAssemble(*vector);
    TrilinosSparseSpaceType::GlobalAssemble(*matrix);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalBuildSystemStructure, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    Teuchos::MpiComm<int> tpetra_comm(raw_mpi_comm);

    const int local_size = 2;
    const int first_my_id = r_comm.Rank() * local_size;
    const int system_size = 2 * r_comm.Size();

    std::vector<std::vector<int>> all_equation_ids = {{first_my_id, first_my_id + 1}};

    TrilinosSparseSpaceType::MatrixPointerType pA;
    TrilinosSparseSpaceType::VectorPointerType pb;
    TrilinosSparseSpaceType::VectorPointerType pDx;
    TrilinosSparseSpaceType::VectorPointerType pReactions;

    TrilinosSparseSpaceType::MapPointerType pMap;
    std::vector<typename TrilinosSparseSpaceType::MapType::global_ordinal_type> local_ids(local_size);
    for (int i = 0; i < local_size; i++) local_ids[i] = first_my_id + i;
    pMap = Teuchos::rcp(new typename TrilinosSparseSpaceType::MapType(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), local_ids, 0, Teuchos::rcp(&tpetra_comm, false)));

    TrilinosSparseSpaceType::BuildSystemStructure(
        tpetra_comm, local_size, first_my_id, 5, all_equation_ids,
        pA, pb, pDx, pReactions, system_size, pMap
    );

    KRATOS_EXPECT_FALSE(TrilinosSparseSpaceType::IsNull(pA));
    KRATOS_EXPECT_FALSE(TrilinosSparseSpaceType::IsNull(pb));
    KRATOS_EXPECT_FALSE(TrilinosSparseSpaceType::IsNull(pDx));
    KRATOS_EXPECT_FALSE(TrilinosSparseSpaceType::IsNull(pReactions));
    KRATOS_EXPECT_EQ(static_cast<std::size_t>(system_size), TrilinosSparseSpaceType::Size1(*pA));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalBuildSystemStructureRowColumnBlocks, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    Teuchos::MpiComm<int> tpetra_comm(raw_mpi_comm);

    const int local_size = 2;
    const int first_my_id = r_comm.Rank() * local_size;
    const int system_size = 2 * r_comm.Size();

    // Rectangular row/column blocks: one diagonal and one cross entry per rank
    std::vector<std::vector<int>> all_row_equation_ids = {{first_my_id}, {first_my_id}};
    std::vector<std::vector<int>> all_col_equation_ids = {{first_my_id}, {first_my_id + 1}};

    TrilinosSparseSpaceType::MatrixPointerType pA;
    TrilinosSparseSpaceType::VectorPointerType pb;
    TrilinosSparseSpaceType::VectorPointerType pDx;
    TrilinosSparseSpaceType::VectorPointerType pReactions;

    TrilinosSparseSpaceType::MapPointerType pMap;
    std::vector<typename TrilinosSparseSpaceType::MapType::global_ordinal_type> local_ids(local_size);
    for (int i = 0; i < local_size; i++) local_ids[i] = first_my_id + i;
    pMap = Teuchos::rcp(new typename TrilinosSparseSpaceType::MapType(
        Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), local_ids, 0,
        Teuchos::rcp(&tpetra_comm, false)));

    TrilinosSparseSpaceType::BuildSystemStructure(
        tpetra_comm, local_size, first_my_id, 5,
        all_row_equation_ids, all_col_equation_ids,
        pA, pb, pDx, pReactions, system_size, pMap
    );

    KRATOS_EXPECT_FALSE(TrilinosSparseSpaceType::IsNull(pA));
    KRATOS_EXPECT_FALSE(TrilinosSparseSpaceType::IsNull(pb));
    KRATOS_EXPECT_FALSE(TrilinosSparseSpaceType::IsNull(pDx));
    KRATOS_EXPECT_FALSE(TrilinosSparseSpaceType::IsNull(pReactions));
    KRATOS_EXPECT_EQ(static_cast<std::size_t>(system_size), TrilinosSparseSpaceType::Size1(*pA));
    // Each rank contributes 2 non-zeros in row first_my_id (cols first_my_id, first_my_id+1)
    KRATOS_EXPECT_EQ(2 * r_comm.Size(), static_cast<int>(pA->getGlobalNumEntries()));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalBuildConstraintsStructure, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    Teuchos::MpiComm<int> tpetra_comm(raw_mpi_comm);

    const int local_size = 2;
    const int first_my_id = r_comm.Rank() * local_size;

    std::vector<std::vector<int>> all_slave_ids = {{first_my_id}};
    std::vector<std::vector<int>> all_master_ids = {{first_my_id + 1}};

    TrilinosSparseSpaceType::MatrixPointerType pT;
    TrilinosSparseSpaceType::VectorPointerType pConstantVector;

    TrilinosSparseSpaceType::MapPointerType pMap;
    std::vector<typename TrilinosSparseSpaceType::MapType::global_ordinal_type> local_ids(local_size);
    for (int i = 0; i < local_size; i++) local_ids[i] = first_my_id + i;
    pMap = Teuchos::rcp(new typename TrilinosSparseSpaceType::MapType(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), local_ids, 0, Teuchos::rcp(&tpetra_comm, false)));

    TrilinosSparseSpaceType::BuildConstraintsStructure(
        tpetra_comm, local_size, first_my_id, 5, all_slave_ids, all_master_ids,
        pT, pConstantVector, pMap
    );

    KRATOS_EXPECT_FALSE(TrilinosSparseSpaceType::IsNull(pT));
    KRATOS_EXPECT_FALSE(TrilinosSparseSpaceType::IsNull(pConstantVector));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalGetCommunicatorMatrix, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    auto matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size);
    auto& r_trilinos_comm = TrilinosSparseSpaceType::GetCommunicator(*matrix);
    KRATOS_EXPECT_EQ(r_comm.Rank(), r_trilinos_comm.getRank());
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalCreateVectorCopy, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    auto vector = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size);
    auto p_copy = TrilinosSparseSpaceType::CreateVectorCopy(*vector);
    auto local_vector = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalVector(size);
    TrilinosCPPTestExperimentalUtilities::CheckSparseVectorFromLocalVector(*p_copy, local_vector);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalCreateEmptyPointerWithCommPtr, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    TrilinosSparseSpaceType::CommunicatorPointerType tpetra_comm =
        Teuchos::rcp(new TrilinosSparseSpaceType::CommunicatorType(raw_mpi_comm));

    // With communicator, size-0 real objects are created (not null)
    auto pA = TrilinosSparseSpaceType::CreateEmptyMatrixPointer(tpetra_comm);
    KRATOS_EXPECT_FALSE(TrilinosSparseSpaceType::IsNull(pA));
    KRATOS_EXPECT_EQ(0, static_cast<int>(TrilinosSparseSpaceType::Size1(*pA)));

    auto pv = TrilinosSparseSpaceType::CreateEmptyVectorPointer(tpetra_comm);
    KRATOS_EXPECT_FALSE(TrilinosSparseSpaceType::IsNull(pv));
    KRATOS_EXPECT_EQ(0, static_cast<int>(TrilinosSparseSpaceType::Size(*pv)));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalCreateEmptyPointerWithCommRef, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    TrilinosSparseSpaceType::CommunicatorType tpetra_comm(raw_mpi_comm);

    // With communicator ref, delegates to CommunicatorPointerType version — size-0 real objects
    auto pA = TrilinosSparseSpaceType::CreateEmptyMatrixPointer(tpetra_comm);
    KRATOS_EXPECT_FALSE(TrilinosSparseSpaceType::IsNull(pA));
    KRATOS_EXPECT_EQ(0, static_cast<int>(TrilinosSparseSpaceType::Size1(*pA)));

    auto pv = TrilinosSparseSpaceType::CreateEmptyVectorPointer(tpetra_comm);
    KRATOS_EXPECT_FALSE(TrilinosSparseSpaceType::IsNull(pv));
    KRATOS_EXPECT_EQ(0, static_cast<int>(TrilinosSparseSpaceType::Size(*pv)));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalCopyVector, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    auto vector_x = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size);
    // vector_y starts with offset=1 values, should be overwritten by Copy
    auto vector_y = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size, 1.0);
    TrilinosSparseSpaceType::Copy(*vector_x, *vector_y);
    auto local_vector = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalVector(size);
    TrilinosCPPTestExperimentalUtilities::CheckSparseVectorFromLocalVector(*vector_y, local_vector);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalCopyMatrix, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    auto matrix_x = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size);
    // matrix_y starts with offset=10 values; after Copy should match matrix_x
    auto matrix_y = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 10.0);
    TrilinosSparseSpaceType::Copy(*matrix_x, *matrix_y);
    auto local_matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalMatrix(size);
    TrilinosCPPTestExperimentalUtilities::CheckSparseMatrixFromLocalMatrix(*matrix_y, local_matrix);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalSetValueVector, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    auto vector = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size);

    // Each rank sets its first owned global index to 42.0
    const int global_idx = r_comm.Rank() * 2;
    TrilinosSparseSpaceType::SetValue(*vector, global_idx, 42.0);

    auto map = vector->getMap();
    const auto local_idx = map->getLocalElement(static_cast<TrilinosSparseSpaceType::GO>(global_idx));
    if (local_idx != Tpetra::Details::OrdinalTraits<TrilinosSparseSpaceType::LO>::invalid()) {
        auto localView = vector->getLocalViewHost(Tpetra::Access::ReadOnly);
        KRATOS_EXPECT_DOUBLE_EQ(42.0, static_cast<double>(localView(local_idx, 0)));
    }
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalSetGlobalVecMethods, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    auto vector = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size, 0.0);

    const int first_my_gid  = r_comm.Rank() * 2;
    const int second_my_gid = first_my_gid + 1;

    // SetGlobalVec: set first owned entry to 42.0 and verify immediately
    TrilinosSparseSpaceType::SetGlobalVec(*vector, first_my_gid, 42.0);
    auto map = vector->getMap();
    auto local_idx = map->getLocalElement(static_cast<TrilinosSparseSpaceType::GO>(first_my_gid));
    if (local_idx != Tpetra::Details::OrdinalTraits<TrilinosSparseSpaceType::LO>::invalid()) {
        auto localView = vector->getLocalViewHost(Tpetra::Access::ReadOnly);
        KRATOS_EXPECT_DOUBLE_EQ(42.0, static_cast<double>(localView(local_idx, 0)));
    }

    // SetGlobalVecNoAssemble: set second owned entry to 21.0 then assemble explicitly
    TrilinosSparseSpaceType::SetGlobalVecNoAssemble(*vector, second_my_gid, 21.0);
    TrilinosSparseSpaceType::GlobalAssemble(*vector);
    local_idx = map->getLocalElement(static_cast<TrilinosSparseSpaceType::GO>(second_my_gid));
    if (local_idx != Tpetra::Details::OrdinalTraits<TrilinosSparseSpaceType::LO>::invalid()) {
        auto localView = vector->getLocalViewHost(Tpetra::Access::ReadOnly);
        KRATOS_EXPECT_DOUBLE_EQ(21.0, static_cast<double>(localView(local_idx, 0)));
    }
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalSetLocalVecMethods, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    auto vector = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size, 0.0);

    // SetLocalVec: set local row 0 to 11.0 and verify
    TrilinosSparseSpaceType::SetLocalVec(*vector, 0, 11.0);
    {
        auto localView = vector->getLocalViewHost(Tpetra::Access::ReadOnly);
        KRATOS_EXPECT_DOUBLE_EQ(11.0, static_cast<double>(localView(0, 0)));
    }

    // SetLocalVecNoAssemble: set local row 1 to 7.0 and assemble explicitly
    TrilinosSparseSpaceType::SetLocalVecNoAssemble(*vector, 1, 7.0);
    TrilinosSparseSpaceType::GlobalAssemble(*vector);
    {
        auto localView = vector->getLocalViewHost(Tpetra::Access::ReadOnly);
        KRATOS_EXPECT_DOUBLE_EQ(7.0, static_cast<double>(localView(1, 0)));
    }
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalSetGlobalMatMethods, KratosTrilinosApplicationMPITestSuite)
{
    using GO = TrilinosSparseSpaceType::GO;
    using LO = TrilinosSparseSpaceType::LO;

    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    const int first_my_gid = r_comm.Rank() * 2;

    // Build a diagonal sparse matrix (only diagonal entries in sparsity pattern)
    auto matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 0.0, false);

    // SetGlobalMat: set diagonal entry at global row=first_my_gid to 5.0
    TrilinosSparseSpaceType::SetGlobalMat(*matrix, first_my_gid, first_my_gid, 5.0);
    {
        auto row_map = matrix->getRowMap();
        LO local_row = row_map->getLocalElement(static_cast<GO>(first_my_gid));
        typename TrilinosSparseSpaceType::MatrixType::local_inds_host_view_type local_cols;
        typename TrilinosSparseSpaceType::MatrixType::values_host_view_type vals;
        matrix->getLocalRowView(local_row, local_cols, vals);
        bool found = false;
        for (std::size_t k = 0; k < local_cols.extent(0); ++k) {
            if (matrix->getColMap()->getGlobalElement(local_cols(k)) == static_cast<GO>(first_my_gid)) {
                KRATOS_EXPECT_DOUBLE_EQ(5.0, static_cast<double>(vals(k)));
                found = true;
            }
        }
        KRATOS_EXPECT_TRUE(found);
    }

    // SetGlobalMatNoAssemble: set diagonal entry at global row=first_my_gid+1 to 3.0
    TrilinosSparseSpaceType::SetGlobalMatNoAssemble(*matrix, first_my_gid + 1, first_my_gid + 1, 3.0);
    TrilinosSparseSpaceType::GlobalAssemble(*matrix);
    {
        auto row_map = matrix->getRowMap();
        LO local_row = row_map->getLocalElement(static_cast<GO>(first_my_gid + 1));
        typename TrilinosSparseSpaceType::MatrixType::local_inds_host_view_type local_cols;
        typename TrilinosSparseSpaceType::MatrixType::values_host_view_type vals;
        matrix->getLocalRowView(local_row, local_cols, vals);
        bool found = false;
        for (std::size_t k = 0; k < local_cols.extent(0); ++k) {
            if (matrix->getColMap()->getGlobalElement(local_cols(k)) == static_cast<GO>(first_my_gid + 1)) {
                KRATOS_EXPECT_DOUBLE_EQ(3.0, static_cast<double>(vals(k)));
                found = true;
            }
        }
        KRATOS_EXPECT_TRUE(found);
    }
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalSetLocalMatMethods, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();

    // Build a diagonal sparse matrix (only diagonal entries in sparsity pattern)
    auto matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 0.0, false);

    // SetLocalMat: set local row 0, local col 0 (diagonal) to 9.0
    TrilinosSparseSpaceType::SetLocalMat(*matrix, 0, 0, 9.0);
    {
        typename TrilinosSparseSpaceType::MatrixType::local_inds_host_view_type local_cols;
        typename TrilinosSparseSpaceType::MatrixType::values_host_view_type vals;
        matrix->getLocalRowView(0, local_cols, vals);
        bool found = false;
        for (std::size_t k = 0; k < local_cols.extent(0); ++k) {
            if (local_cols(k) == 0) {
                KRATOS_EXPECT_DOUBLE_EQ(9.0, static_cast<double>(vals(k)));
                found = true;
            }
        }
        KRATOS_EXPECT_TRUE(found);
    }

    // SetLocalMatNoAssemble: set local row 1, local col 1 (diagonal) to 4.0
    TrilinosSparseSpaceType::SetLocalMatNoAssemble(*matrix, 1, 1, 4.0);
    TrilinosSparseSpaceType::GlobalAssemble(*matrix);
    {
        typename TrilinosSparseSpaceType::MatrixType::local_inds_host_view_type local_cols;
        typename TrilinosSparseSpaceType::MatrixType::values_host_view_type vals;
        matrix->getLocalRowView(1, local_cols, vals);
        bool found = false;
        for (std::size_t k = 0; k < local_cols.extent(0); ++k) {
            if (local_cols(k) == 1) {
                KRATOS_EXPECT_DOUBLE_EQ(4.0, static_cast<double>(vals(k)));
                found = true;
            }
        }
        KRATOS_EXPECT_TRUE(found);
    }
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalGetValue, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    auto vector = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size);
    const int first_my_id = r_comm.Rank() * 2;
    KRATOS_EXPECT_DOUBLE_EQ(static_cast<double>(first_my_id),
        TrilinosSparseSpaceType::GetValue(*vector, first_my_id));
    KRATOS_EXPECT_DOUBLE_EQ(static_cast<double>(first_my_id + 1),
        TrilinosSparseSpaceType::GetValue(*vector, first_my_id + 1));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalGatherValues, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    auto vector = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size);
    // Gather the first two global values on every process
    std::vector<int> indices = {0, 1};
    std::vector<double> values(2, 0.0);
    TrilinosSparseSpaceType::GatherValues(*vector, indices, values.data());
    KRATOS_EXPECT_DOUBLE_EQ(0.0, values[0]);
    KRATOS_EXPECT_DOUBLE_EQ(1.0, values[1]);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalClear, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    TrilinosSparseSpaceType::CommunicatorPointerType tpetra_comm =
        Teuchos::rcp(new TrilinosSparseSpaceType::CommunicatorType(raw_mpi_comm));
    const int size = 2 * r_comm.Size();

    auto pA = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size);
    KRATOS_EXPECT_FALSE(TrilinosSparseSpaceType::IsNull(pA));
    TrilinosSparseSpaceType::Clear(pA);
    KRATOS_EXPECT_FALSE(TrilinosSparseSpaceType::IsNull(pA)); // still valid after swap
    KRATOS_EXPECT_EQ(0, static_cast<int>(TrilinosSparseSpaceType::Size1(*pA)));

    auto pv = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size);
    KRATOS_EXPECT_FALSE(TrilinosSparseSpaceType::IsNull(pv));
    TrilinosSparseSpaceType::Clear(pv);
    KRATOS_EXPECT_FALSE(TrilinosSparseSpaceType::IsNull(pv)); // still valid after swap
    KRATOS_EXPECT_EQ(0, static_cast<int>(TrilinosSparseSpaceType::Size(*pv)));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalFastestDirectSolverList, KratosTrilinosApplicationMPITestSuite)
{
    const auto solvers = TrilinosSparseSpaceType::FastestDirectSolverList();
    KRATOS_EXPECT_GT(solvers.size(), std::size_t(0));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalGetDiagonalNorm, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 12;
    auto matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 1.0);
    // Diagonal values are 1, 2, ..., 12 → norm = sqrt(1² + 2² + ... + 12²)
    double ref = 0.0;
    for (int i = 1; i <= size; ++i) ref += static_cast<double>(i) * static_cast<double>(i);
    ref = std::sqrt(ref);
    KRATOS_EXPECT_NEAR(TrilinosSparseSpaceType::GetDiagonalNorm(*matrix), ref, 1.0e-6);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalGetMaxDiagonal, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 12;
    auto matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(12.0, TrilinosSparseSpaceType::GetMaxDiagonal(*matrix));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalAssembleLHS, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    const int first_my_id = r_comm.Rank() * 2;
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    TrilinosSparseSpaceType::CommunicatorPointerType tpetra_comm =
        Teuchos::rcp(new TrilinosSparseSpaceType::CommunicatorType(raw_mpi_comm));

    using LO = TrilinosSparseSpaceType::LO;
    using GO = TrilinosSparseSpaceType::GO;

    // Build FECrsGraph with full 2x2 dense block per rank
    auto map = Teuchos::rcp(new TrilinosSparseSpaceType::MapType(size, 0, tpetra_comm));
    auto graph = Teuchos::rcp(new TrilinosSparseSpaceType::GraphType(map, map, 2));
    std::vector<GO> col_gids = {static_cast<GO>(first_my_id), static_cast<GO>(first_my_id + 1)};
    for (int gid = first_my_id; gid < first_my_id + 2; ++gid) {
        const GO global_gid = static_cast<GO>(gid);
        graph->insertGlobalIndices(global_gid, Teuchos::ArrayView<const GO>(col_gids.data(), 2));
    }
    graph->fillComplete();

    TrilinosSparseMatrixType A(graph);
    A.beginAssembly();

    // Assemble 2x2 local diagonal-like matrix
    Matrix local_mat(2, 2, 0.0);
    local_mat(0, 0) = 5.0;
    local_mat(1, 1) = 3.0;
    std::vector<std::size_t> eq_ids = {static_cast<std::size_t>(first_my_id),
                                        static_cast<std::size_t>(first_my_id + 1)};
    TrilinosSparseSpaceType::AssembleLHS(A, local_mat, eq_ids);
    A.endAssembly();
    if (A.isFillActive()) A.fillComplete();

    // Check diagonal entries
    const LO local_row_0 = map->getLocalElement(static_cast<GO>(first_my_id));
    if (local_row_0 != Tpetra::Details::OrdinalTraits<LO>::invalid()) {
        auto localMatrix = A.getLocalMatrixHost();
        auto localRow0 = localMatrix.row(local_row_0);
        // Find the diagonal (column == first_my_id)
        for (int k = 0; k < localRow0.length; ++k) {
            const GO col_gid = A.getColMap()->getGlobalElement(localRow0.colidx(k));
            if (col_gid == static_cast<GO>(first_my_id)) {
                KRATOS_EXPECT_DOUBLE_EQ(5.0, static_cast<double>(localRow0.value(k)));
            }
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalAssembleRHS, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    const int first_my_id = r_comm.Rank() * 2;
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    TrilinosSparseSpaceType::CommunicatorPointerType tpetra_comm =
        Teuchos::rcp(new TrilinosSparseSpaceType::CommunicatorType(raw_mpi_comm));

    auto map = Teuchos::rcp(new TrilinosSparseSpaceType::MapType(size, 0, tpetra_comm));
    auto p_b = TrilinosSparseSpaceType::CreateVector(map);
    auto& b = *p_b;

    auto p_fe_b = dynamic_cast<TrilinosVectorType*>(&b);
    if (p_fe_b) p_fe_b->beginAssembly();

    Vector local_vec(2);
    local_vec[0] = 7.0;
    local_vec[1] = 3.0;
    std::vector<std::size_t> eq_ids = {static_cast<std::size_t>(first_my_id),
                                        static_cast<std::size_t>(first_my_id + 1)};
    TrilinosSparseSpaceType::AssembleRHS(b, local_vec, eq_ids);

    if (p_fe_b) p_fe_b->endAssembly();

    using GO = TrilinosSparseSpaceType::GO;
    using LO = TrilinosSparseSpaceType::LO;
    auto localView = b.getLocalViewHost(Tpetra::Access::ReadOnly);
    const LO local_0 = map->getLocalElement(static_cast<GO>(first_my_id));
    const LO local_1 = map->getLocalElement(static_cast<GO>(first_my_id + 1));
    if (local_0 != Tpetra::Details::OrdinalTraits<LO>::invalid()) {
        KRATOS_EXPECT_DOUBLE_EQ(7.0, static_cast<double>(localView(local_0, 0)));
    }
    if (local_1 != Tpetra::Details::OrdinalTraits<LO>::invalid()) {
        KRATOS_EXPECT_DOUBLE_EQ(3.0, static_cast<double>(localView(local_1, 0)));
    }
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalCreateDofUpdater, KratosTrilinosApplicationMPITestSuite)
{
    auto dof_updater = TrilinosSparseSpaceType::CreateDofUpdater();
    KRATOS_EXPECT_NE(nullptr, dof_updater.get());
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalGetColumnThrows, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    auto matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size);
    auto p_col = TrilinosSparseSpaceType::CreateVector(matrix->getDomainMap());
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        TrilinosSparseSpaceType::GetColumn(0, *matrix, *p_col),
        "GetColumn method is not currently implemented"
    );
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalResizeThrows, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    auto matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size);
    auto vector = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        TrilinosSparseSpaceType::Resize(*matrix, size, size),
        "Resize is not defined for Trilinos Sparse Matrix"
    );
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        TrilinosSparseSpaceType::Resize(*vector, size),
        "Resize is not defined for a reference to Trilinos Vector"
    );
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalCreateEmptyMapPointer, KratosTrilinosApplicationMPITestSuite)
{
    // CreateEmptyMapPointer always returns Teuchos::null (a null RCP)
    auto p_map = TrilinosSparseSpaceType::CreateEmptyMapPointer();
    KRATOS_EXPECT_TRUE(p_map == Teuchos::null);
    KRATOS_EXPECT_TRUE(TrilinosSparseSpaceType::IsNull(p_map));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalGetMinDiagonal, KratosTrilinosApplicationMPITestSuite)
{
    // Diagonal entries of the dummy matrix (scale=1) are 1, 2, ..., size → min = 1.0
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 12;
    auto matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 1.0);
    // min diagonal is the smallest entry = 1.0 (global row 0 has diagonal = scale * 1)
    const double diag_min = TrilinosSparseSpaceType::GetMinDiagonal(*matrix);
    KRATOS_EXPECT_NEAR(diag_min, 1.0, 1.0e-10);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalGetAveragevalueDiagonal, KratosTrilinosApplicationMPITestSuite)
{
    // Diagonal entries are 1, 2, ..., 12 → average = (1+12)/2 = 6.5
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 12;
    auto matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 1.0);
    const double expected_avg = 0.5 * (1.0 + 12.0);
    KRATOS_EXPECT_NEAR(TrilinosSparseSpaceType::GetAveragevalueDiagonal(*matrix), expected_avg, 1.0e-10);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalManualFinalize, KratosTrilinosApplicationMPITestSuite)
{
    // ManualFinalize calls endAssembly() on a matrix that's open for FE assembly
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    auto matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size);
    // Open the matrix for FE assembly (this is valid after the matrix is in closed state)
    matrix->beginAssembly();
    KRATOS_EXPECT_TRUE(matrix->isFillActive());
    // ManualFinalize must close it without throwing
    TrilinosSparseSpaceType::ManualFinalize(*matrix);
    KRATOS_EXPECT_FALSE(matrix->isFillActive());
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalResizeVectorPointer, KratosTrilinosApplicationMPITestSuite)
{
    // Resize(VectorPointerType&, n) replaces the pointed-to vector with a new empty one of global size n
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int original_size = 2 * r_comm.Size();
    auto p_vec = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, original_size);
    KRATOS_EXPECT_EQ(TrilinosSparseSpaceType::Size(*p_vec), static_cast<std::size_t>(original_size));
    // Resize to 0 — pointer must be updated in-place (pass-by-ref)
    TrilinosSparseSpaceType::Resize(p_vec, std::size_t(0));
    KRATOS_EXPECT_EQ(TrilinosSparseSpaceType::Size(*p_vec), std::size_t(0));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalGetOrCreateMap, KratosTrilinosApplicationMPITestSuite)
{
    // GetOrCreateMap creates a Tpetra::Map with LocalSize entries starting at FirstMyId
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    TrilinosSparseSpaceType::CommunicatorType tpetra_comm(raw_mpi_comm);

    const int rank = r_comm.Rank();
    const int local_size = 2;
    const int first_my_id = rank * local_size;
    auto p_map = TrilinosSparseSpaceType::GetOrCreateMap(tpetra_comm, local_size, first_my_id);

    KRATOS_EXPECT_NE(p_map, Teuchos::null);
    KRATOS_EXPECT_EQ(static_cast<int>(p_map->getNodeNumElements()), local_size);
    // Check that global IDs for this rank start at first_my_id
    auto node_elements = p_map->getNodeElementList();
    KRATOS_EXPECT_EQ(static_cast<int>(node_elements[0]), first_my_id);
    KRATOS_EXPECT_EQ(static_cast<int>(node_elements[1]), first_my_id + 1);
}

} // namespace Kratos::Testing