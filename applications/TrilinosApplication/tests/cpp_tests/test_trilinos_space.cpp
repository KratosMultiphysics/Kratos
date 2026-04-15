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
#include <Teuchos_RCP.hpp>

// Project includes
#include "tests/cpp_tests/trilinos_fast_suite.h"
#include "trilinos_space.h"
#include "containers/model.h"
#include "mpi/includes/mpi_data_communicator.h"
#include "trilinos_cpp_test_utilities.h"
#include "utilities/math_utils.h"

namespace Kratos::Testing
{

/// Basic definitions
using TrilinosSparseMatrix = Epetra_FECrsMatrix;
using TrilinosSparseVector = Epetra_FEVector;
using TrilinosSparseSpaceType = TrilinosSpace<TrilinosSparseMatrix, TrilinosSparseVector>;
using TrilinosLocalSpaceType = UblasSpace<double, Matrix, Vector>;

using TrilinosSparseMatrixType = TrilinosSparseSpaceType::MatrixType;
using TrilinosVectorType = TrilinosSparseSpaceType::VectorType;

using TrilinosLocalMatrixType = TrilinosLocalSpaceType::MatrixType;
using TrilinosLocalVectorType = TrilinosLocalSpaceType::VectorType;

KRATOS_TEST_CASE_IN_SUITE(TrilinosSizeVector, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto vector = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size);

    // Check
    KRATOS_EXPECT_EQ(static_cast<std::size_t>(size), TrilinosSparseSpaceType::Size(vector));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosSizeMatrix, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto matrix = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size);

    // Check
    KRATOS_EXPECT_EQ(static_cast<std::size_t>(size), TrilinosSparseSpaceType::Size1(matrix));
    KRATOS_EXPECT_EQ(static_cast<std::size_t>(size), TrilinosSparseSpaceType::Size2(matrix));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosDotProduct, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto vector1 = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size);
    auto vector2 = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size);
    double ref = 0.0;
    for (int i = 0; i < size; ++i) {
        ref += std::pow(static_cast<double>(i), 2);
    }

    // Check
    KRATOS_EXPECT_DOUBLE_EQ(ref, TrilinosSparseSpaceType::Dot(vector1, vector2));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosMaxMin, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto vector = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size);

    // Check
    KRATOS_EXPECT_DOUBLE_EQ(0.0, TrilinosSparseSpaceType::Min(vector));
    KRATOS_EXPECT_DOUBLE_EQ(static_cast<double>(size - 1), TrilinosSparseSpaceType::Max(vector));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosTwoNormVector, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto vector = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size);
    double ref = 0.0;
    for (int i = 0; i < size; ++i) {
        ref += std::pow(static_cast<double>(i), 2);
    }
    ref = std::sqrt(ref);

    // Check
    KRATOS_EXPECT_DOUBLE_EQ(ref, TrilinosSparseSpaceType::TwoNorm(vector));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosTwoNormMatrix1, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 2 * r_comm.Size();
    auto matrix = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size);
    double ref = 0.0;
    for (int i = 0; i < size; ++i) {
        ref += std::pow(static_cast<double>(i), 2);
    }
    ref = std::sqrt(ref);

    // Check
    KRATOS_EXPECT_DOUBLE_EQ(ref, TrilinosSparseSpaceType::TwoNorm(matrix));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosTwoNormMatrix2, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 2 * r_comm.Size();
    auto matrix = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 0.0, true);
    auto local_matrix = TrilinosCPPTestUtilities::GenerateDummyLocalMatrix(size, 0.0, true);

    // Check
    KRATOS_EXPECT_DOUBLE_EQ(TrilinosLocalSpaceType::TwoNorm(local_matrix), TrilinosSparseSpaceType::TwoNorm(matrix));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosMultMatrixVector, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 2 * r_comm.Size();
    auto matrix = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 0.0, true);
    auto local_matrix = TrilinosCPPTestUtilities::GenerateDummyLocalMatrix(size, 0.0, true);
    auto vector = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size, 0.0);
    auto local_vector = TrilinosCPPTestUtilities::GenerateDummyLocalVector(size, 0.0);

    // Epetra communicator
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    TrilinosSparseSpaceType::CommunicatorType epetra_comm(raw_mpi_comm);

    // Create a map
    TrilinosSparseSpaceType::MapType Map(size,0,epetra_comm);

    // Create an Epetra_Vector
    TrilinosVectorType mult(Map);

    // Solution
    TrilinosSparseSpaceType::Mult(matrix, vector, mult);

    // Check
    const TrilinosLocalVectorType multiply_reference = prod(local_matrix, local_vector);
    TrilinosCPPTestUtilities::CheckSparseVectorFromLocalVector(mult, multiply_reference);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosMultMatrixMatrix, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 2 * r_comm.Size();
    auto matrix_1 = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 0.0, true);
    auto local_matrix_1 = TrilinosCPPTestUtilities::GenerateDummyLocalMatrix(size, 0.0, true);
    auto matrix_2 = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 1.0, true);
    auto local_matrix_2 = TrilinosCPPTestUtilities::GenerateDummyLocalMatrix(size, 1.0, true);

    // Epetra communicator
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    TrilinosSparseSpaceType::CommunicatorType epetra_comm(raw_mpi_comm);

    // Create a map
    TrilinosSparseSpaceType::MapType Map(size,0,epetra_comm);

    // Create an Epetra_Matrix
    std::vector<int> NumNz;
    TrilinosSparseMatrixType mult(Copy, Map, NumNz.data());

    // Solution
    TrilinosSparseSpaceType::Mult(matrix_1, matrix_2, mult);

    // Check
    const TrilinosLocalMatrixType multiply_reference = prod(local_matrix_1, local_matrix_2);
    TrilinosCPPTestUtilities::CheckSparseMatrixFromLocalMatrix(mult, multiply_reference);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosTransposeMultMatrixVector, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 2 * r_comm.Size();
    auto matrix = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 0.0, true);
    auto local_matrix = TrilinosCPPTestUtilities::GenerateDummyLocalMatrix(size, 0.0, true);
    auto vector = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size, 0.0);
    auto local_vector = TrilinosCPPTestUtilities::GenerateDummyLocalVector(size, 0.0);

    // Epetra communicator
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    TrilinosSparseSpaceType::CommunicatorType epetra_comm(raw_mpi_comm);

    // Create a map
    TrilinosSparseSpaceType::MapType Map(size,0,epetra_comm);

    // Create an Epetra_Vector
    TrilinosVectorType mult(Map);

    // Solution
    TrilinosSparseSpaceType::TransposeMult(matrix, vector, mult);

    // Check
    const TrilinosLocalVectorType multiply_reference = prod(trans(local_matrix), local_vector);
    TrilinosCPPTestUtilities::CheckSparseVectorFromLocalVector(mult, multiply_reference);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosTransposeMultMatrixMatrix, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 2 * r_comm.Size();
    auto matrix_1 = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 0.0, true);
    auto local_matrix_1 = TrilinosCPPTestUtilities::GenerateDummyLocalMatrix(size, 0.0, true);
    auto matrix_2 = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 1.0, true);
    auto local_matrix_2 = TrilinosCPPTestUtilities::GenerateDummyLocalMatrix(size, 1.0, true);

    // Epetra communicator
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    TrilinosSparseSpaceType::CommunicatorType epetra_comm(raw_mpi_comm);

    // Create a map
    TrilinosSparseSpaceType::MapType Map(size,0,epetra_comm);

    // Create an Epetra_Matrix
    std::vector<int> NumNz;
    TrilinosSparseMatrixType mult(Copy, Map, NumNz.data());

    // Solution
    TrilinosSparseSpaceType::TransposeMult(matrix_1, matrix_2, mult, {true, false});

    // Check
    const TrilinosLocalMatrixType multiply_reference = prod(trans(local_matrix_1), local_matrix_2);
    TrilinosCPPTestUtilities::CheckSparseMatrixFromLocalMatrix(mult, multiply_reference);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosBtDBProductOperation, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 2 * r_comm.Size();
    auto matrix_1 = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 0.0, true);
    auto local_matrix_1 = TrilinosCPPTestUtilities::GenerateDummyLocalMatrix(size, 0.0, true);
    auto matrix_2 = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 1.0, true);
    auto local_matrix_2 = TrilinosCPPTestUtilities::GenerateDummyLocalMatrix(size, 1.0, true);

    // Epetra communicator
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    TrilinosSparseSpaceType::CommunicatorType epetra_comm(raw_mpi_comm);

    // Create a map
    TrilinosSparseSpaceType::MapType Map(size,0,epetra_comm);

    // Create an Epetra_Matrix
    std::vector<int> NumNz;
    TrilinosSparseMatrixType mult(Copy, Map, NumNz.data());

    // Solution
    TrilinosSparseSpaceType::BtDBProductOperation(mult, matrix_1, matrix_2);

    // Check
    TrilinosLocalMatrixType multiply_reference;
    MathUtils<double>::BtDBProductOperation(multiply_reference, local_matrix_1, local_matrix_2);
    TrilinosCPPTestUtilities::CheckSparseMatrixFromLocalMatrix(mult, multiply_reference);

    // Non zero matrix
    auto second_mult = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 0.0, true);

    // Solution
    TrilinosSparseSpaceType::BtDBProductOperation(second_mult, matrix_1, matrix_2);
    TrilinosCPPTestUtilities::CheckSparseMatrixFromLocalMatrix(second_mult, multiply_reference);
}

// Error related to Trilinos issue: https://github.com/trilinos/Trilinos/issues/9252
// KRATOS_TEST_CASE_IN_SUITE(TrilinosBtDBProductOperationRealCase, KratosTrilinosApplicationMPITestSuite)
// {
//     // The data communicator
//     const auto& r_comm = Testing::GetDefaultDataCommunicator();

//     // The dummy matrix
//     const int size = 6;

//     // Generate A matrix
//     std::vector<int> row_indexes = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5};
//     std::vector<int> column_indexes = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 2, 3, 4, 5, 2, 3, 4, 5};
//     std::vector<double> values = {2069000000.0, 0.0, -2069000000.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2069000000.0, 0.0, 4138000000.0, 0.0, -2069000000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2069000000.0, 0.0, 2069000000.0, 0.0, 0.0, 0.0, 0.0, 0.0};
//     auto A = TrilinosCPPTestUtilities::GenerateSparseMatrix(r_comm, size, row_indexes, column_indexes, values);
//     const TrilinosSparseSpaceType::MapType* p_map = &(A.RowMap());

//     // Generate T matrix
//     row_indexes = {0,1,2,3,4,4,5};
//     column_indexes = {0,1,2,3,2,4,5};
//     values = {1,1,1,1,1,0,1};
//     auto T = TrilinosCPPTestUtilities::GenerateSparseMatrix(r_comm, size, row_indexes, column_indexes, values, p_map);

//     /* Intermediate multiplication */

//     // Create an Epetra_Matrix
//     TrilinosSparseMatrixType aux(::Copy, A.Graph());

//     // First multiplication
//     TrilinosSparseSpaceType::TransposeMult(T, A, aux, {true, false}, true, true);

//     // Values to check
//     row_indexes = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5};
//     column_indexes = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 2, 3, 4, 5, 2, 3, 4, 5};
//     values = {2069000000.0, 0.0, -2069000000.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2069000000.0, 0.0, 2069000000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

//     // Check
//     TrilinosCPPTestUtilities::CheckSparseMatrix(aux, row_indexes, column_indexes, values);

//     // Compute T' A T
//     const TrilinosSparseMatrixType copy_A(A);
//     TrilinosSparseSpaceType::BtDBProductOperation(A, copy_A, T, true, false, true);

//     // Values to check
//     row_indexes = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5};
//     column_indexes = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 2, 3, 4, 5, 2, 3, 4, 5};
//     values = {2069000000.0, 0.0, -2069000000.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2069000000.0, 0.0, 2069000000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

//     // Check
//     TrilinosCPPTestUtilities::CheckSparseMatrix(A, row_indexes, column_indexes, values);
// }

KRATOS_TEST_CASE_IN_SUITE(TrilinosBDBtProductOperation, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 2 * r_comm.Size();
    auto matrix_1 = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 0.0, true);
    auto local_matrix_1 = TrilinosCPPTestUtilities::GenerateDummyLocalMatrix(size, 0.0, true);
    auto matrix_2 = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 1.0, true);
    auto local_matrix_2 = TrilinosCPPTestUtilities::GenerateDummyLocalMatrix(size, 1.0, true);

    // Epetra communicator
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    TrilinosSparseSpaceType::CommunicatorType epetra_comm(raw_mpi_comm);

    // Create a map
    TrilinosSparseSpaceType::MapType Map(size,0,epetra_comm);

    // Create an Epetra_Matrix
    std::vector<int> NumNz;
    TrilinosSparseMatrixType mult(Copy, Map, NumNz.data());

    // Solution
    TrilinosSparseSpaceType::BDBtProductOperation(mult, matrix_1, matrix_2);

    // Check
    TrilinosLocalMatrixType multiply_reference;
    MathUtils<double>::BDBtProductOperation(multiply_reference, local_matrix_1, local_matrix_2);
    TrilinosCPPTestUtilities::CheckSparseMatrixFromLocalMatrix(mult, multiply_reference);

    // Non zero matrix
    auto second_mult = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 0.0, true);

    // Solution
    TrilinosSparseSpaceType::BDBtProductOperation(second_mult, matrix_1, matrix_2);
    TrilinosCPPTestUtilities::CheckSparseMatrixFromLocalMatrix(second_mult, multiply_reference);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosInplaceMult, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto vector = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size);
    auto local_vector = TrilinosCPPTestUtilities::GenerateDummyLocalVector(size);

    // Multiply
    const double mult = 2.0;

    // Solution
    TrilinosSparseSpaceType::InplaceMult(vector, mult);

    // Check
    local_vector *= mult;
    TrilinosCPPTestUtilities::CheckSparseVectorFromLocalVector(vector, local_vector);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosAssign, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto vector_1 = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size);
    auto vector_2 = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size, 1.0);
    auto local_vector_1 = TrilinosCPPTestUtilities::GenerateDummyLocalVector(size);
    auto local_vector_2 = TrilinosCPPTestUtilities::GenerateDummyLocalVector(size, 1.0);

    // Multiply
    const double mult = 2.0;

    // Solution
    TrilinosSparseSpaceType::Assign(vector_1, mult, vector_2);

    // Check
    local_vector_1 = mult * local_vector_2;
    TrilinosCPPTestUtilities::CheckSparseVectorFromLocalVector(vector_1, local_vector_1);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosUnaliasedAdd, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto vector_1 = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size);
    auto vector_2 = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size, 1.0);
    auto local_vector_1 = TrilinosCPPTestUtilities::GenerateDummyLocalVector(size);
    auto local_vector_2 = TrilinosCPPTestUtilities::GenerateDummyLocalVector(size, 1.0);

    // Multiply
    const double mult = 2.0;

    // Solution
    TrilinosSparseSpaceType::UnaliasedAdd(vector_1, mult, vector_2);

    // Check
    local_vector_1 += mult * local_vector_2;
    TrilinosCPPTestUtilities::CheckSparseVectorFromLocalVector(vector_1, local_vector_1);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosScaleAndAdd1, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto vector_1 = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size);
    auto vector_2 = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size, 1.0);
    auto vector_3 = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size, 2.0);
    auto local_vector_1 = TrilinosCPPTestUtilities::GenerateDummyLocalVector(size);
    auto local_vector_2 = TrilinosCPPTestUtilities::GenerateDummyLocalVector(size, 1.0);
    auto local_vector_3 = TrilinosCPPTestUtilities::GenerateDummyLocalVector(size, 2.0);

    // Multiply
    const double mult_1 = 2.0;
    const double mult_2 = 1.5;

    // Solution
    TrilinosSparseSpaceType::ScaleAndAdd(mult_1, vector_2, mult_2, vector_3, vector_1);

    // Check
    local_vector_1 = mult_1 * local_vector_2 + mult_2 * local_vector_3;
    TrilinosCPPTestUtilities::CheckSparseVectorFromLocalVector(vector_1, local_vector_1);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosScaleAndAdd2, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto vector_1 = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size);
    auto vector_2 = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size, 1.0);
    auto local_vector_1 = TrilinosCPPTestUtilities::GenerateDummyLocalVector(size);
    auto local_vector_2 = TrilinosCPPTestUtilities::GenerateDummyLocalVector(size, 1.0);

    // Multiply
    const double mult_1 = 2.0;
    const double mult_2 = 1.5;

    // Solution
    TrilinosSparseSpaceType::ScaleAndAdd(mult_1, vector_2, mult_2, vector_1);

    // Check
    local_vector_1 = mult_1 * local_vector_2 + mult_2 * local_vector_1;
    TrilinosCPPTestUtilities::CheckSparseVectorFromLocalVector(vector_1, local_vector_1);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosSet, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto vector = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size);
    auto local_vector = TrilinosCPPTestUtilities::GenerateDummyLocalVector(size);

    // Set
    const double value = 2.0;

    // Solution
    TrilinosSparseSpaceType::Set(vector, value);

    // Check
    for (int i = 0; i < size; ++i) local_vector[i] = value;
    TrilinosCPPTestUtilities::CheckSparseVectorFromLocalVector(vector, local_vector);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosSetToZeroMatrix, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto matrix = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size);
    TrilinosLocalMatrixType local_matrix = ZeroMatrix(size, size);

    // Solution
    TrilinosSparseSpaceType::SetToZero(matrix);

    // Check
    TrilinosCPPTestUtilities::CheckSparseMatrixFromLocalMatrix(matrix, local_matrix);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosSetToZeroVector, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto vector = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size);
    TrilinosLocalVectorType local_vector = ZeroVector(size);

    // Solution
    TrilinosSparseSpaceType::SetToZero(vector);

    // Check
    TrilinosCPPTestUtilities::CheckSparseVectorFromLocalVector(vector, local_vector);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosCopyMatrixValues, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto matrix_1 = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 0.0, true);
    auto matrix_2 = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 0.0);
    auto local_matrix = TrilinosCPPTestUtilities::GenerateDummyLocalMatrix(size, 0.0);

    // Solution
    TrilinosSparseSpaceType::CopyMatrixValues(matrix_1, matrix_2);

    // Check
    TrilinosCPPTestUtilities::CheckSparseMatrixFromLocalMatrix(matrix_1, local_matrix);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosCombineMatricesGraphs, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto matrix_1 = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 0.0);
    auto matrix_2 = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 0.0, true);
    auto local_matrix = TrilinosCPPTestUtilities::GenerateDummyLocalMatrix(size, 0.0, true);

    // Creating new matrix
    const auto combined_graph = TrilinosSparseSpaceType::CombineMatricesGraphs(matrix_1, matrix_2);
    TrilinosSparseMatrixType copied_matrix(::Copy, combined_graph);

    // Solution
    TrilinosSparseSpaceType::CopyMatrixValues(copied_matrix, matrix_2);

    // Check
    TrilinosCPPTestUtilities::CheckSparseMatrixFromLocalMatrix(copied_matrix, local_matrix);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosCheckAndCorrectZeroDiagonalValues, KratosTrilinosApplicationMPITestSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");
    auto& r_process_info = r_model_part.GetProcessInfo();
    r_process_info.SetValue(BUILD_SCALE_FACTOR, 1.0);

    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 12;
    auto matrix12x12 = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size);

    // Generate Epetra communicator
    KRATOS_ERROR_IF_NOT(r_comm.IsDistributed()) << "Only distributed DataCommunicators can be used!" << std::endl;
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    TrilinosSparseSpaceType::CommunicatorType epetra_comm(raw_mpi_comm);

    // Dummy vector
    TrilinosSparseSpaceType::MapType map(size, 0, epetra_comm);
    TrilinosVectorType vector12(map);

    // Test the norm of the matrix
    double norm = 0.0;
    norm = TrilinosSparseSpaceType::CheckAndCorrectZeroDiagonalValues(r_process_info, matrix12x12, vector12, SCALING_DIAGONAL::NO_SCALING);
    KRATOS_EXPECT_DOUBLE_EQ(norm, 1.0);
    if (r_comm.Rank() == 0) {
        const auto& raw_values = matrix12x12.ExpertExtractValues();
        KRATOS_EXPECT_DOUBLE_EQ(raw_values[0], 1.0);
    }
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosIsDistributed, KratosTrilinosApplicationMPITestSuite)
{
    KRATOS_EXPECT_TRUE(TrilinosSparseSpaceType::IsDistributed());
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosGetScaleNorm, KratosTrilinosApplicationMPITestSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");
    auto& r_process_info = r_model_part.GetProcessInfo();
    r_process_info.SetValue(BUILD_SCALE_FACTOR, 3.0);

    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 12;
    auto matrix12x12 = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 1.0);

    // Test the norm of the matrix
    double norm = 0.0;
    norm = TrilinosSparseSpaceType::GetScaleNorm(r_process_info, matrix12x12, SCALING_DIAGONAL::NO_SCALING);
    KRATOS_EXPECT_DOUBLE_EQ(norm, 1.0);
    norm = TrilinosSparseSpaceType::GetScaleNorm(r_process_info, matrix12x12, SCALING_DIAGONAL::CONSIDER_PRESCRIBED_DIAGONAL);
    KRATOS_EXPECT_DOUBLE_EQ(norm, 3.0);
    norm = TrilinosSparseSpaceType::GetScaleNorm(r_process_info, matrix12x12, SCALING_DIAGONAL::CONSIDER_NORM_DIAGONAL);
    KRATOS_EXPECT_NEAR(norm, 2.124591464, 1.0e-6);
    norm = TrilinosSparseSpaceType::GetScaleNorm(r_process_info, matrix12x12, SCALING_DIAGONAL::CONSIDER_MAX_DIAGONAL);
    KRATOS_EXPECT_DOUBLE_EQ(norm, 12.0);
    norm = TrilinosSparseSpaceType::GetAveragevalueDiagonal(matrix12x12);
    KRATOS_EXPECT_DOUBLE_EQ(norm, 6.5);
    norm = TrilinosSparseSpaceType::GetMinDiagonal(matrix12x12);
    KRATOS_EXPECT_DOUBLE_EQ(norm, 1.0);
}


KRATOS_TEST_CASE_IN_SUITE(TrilinosScaleAndAddMatrix, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 12;
    auto matrix_1 = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 10.0, true);
    auto matrix_2 = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 20.0, true);

    // Solution
    TrilinosSparseSpaceType::ScaleAndAdd(2.0, matrix_1, 3.0, matrix_2);

    // Check
    auto local_matrix_1 = TrilinosCPPTestUtilities::GenerateDummyLocalMatrix(size, 10.0, true);
    auto local_matrix_2 = TrilinosCPPTestUtilities::GenerateDummyLocalMatrix(size, 20.0, true);
    TrilinosLocalMatrixType local_reference = 2.0 * local_matrix_1 + 3.0 * local_matrix_2;

    TrilinosCPPTestUtilities::CheckSparseMatrixFromLocalMatrix(matrix_2, local_reference);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosMatrixMarket, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 12;
    auto matrix = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 10.0, true);

    // Write to file
    const std::string file_name = "test_matrix_epetra.mm";
    TrilinosSparseSpaceType::WriteMatrixMarketMatrix(file_name.c_str(), matrix, false);

    // Epetra communicator
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    TrilinosSparseSpaceType::CommunicatorType epetra_comm(raw_mpi_comm);

    // Read from file
    auto read_matrix = TrilinosSparseSpaceType::ReadMatrixMarket(file_name, epetra_comm);

    // Check
    auto local_matrix = TrilinosCPPTestUtilities::GenerateDummyLocalMatrix(size, 10.0, true);
    TrilinosCPPTestUtilities::CheckSparseMatrixFromLocalMatrix(*read_matrix, local_matrix);

    // Clean up
    if (r_comm.Rank() == 0) {
        std::remove(file_name.c_str());
    }
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosMatrixMarketVector, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 12;
    auto vector = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size);

    // Write to file
    const std::string file_name = "test_vector_epetra.mm";
    TrilinosSparseSpaceType::WriteMatrixMarketVector(file_name.c_str(), vector);

    // Epetra communicator
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    TrilinosSparseSpaceType::CommunicatorType epetra_comm(raw_mpi_comm);

    // Read from file
    auto read_vector = TrilinosSparseSpaceType::ReadMatrixMarketVector(file_name, epetra_comm, size);

    // Check
    auto local_vector = TrilinosCPPTestUtilities::GenerateDummyLocalVector(size);
    TrilinosCPPTestUtilities::CheckSparseVectorFromLocalVector(*read_vector, local_vector);

    // Clean up
    if (r_comm.Rank() == 0) {
        std::remove(file_name.c_str());
    }
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosIsDistributedSpace, KratosTrilinosApplicationMPITestSuite)
{
    KRATOS_EXPECT_TRUE(TrilinosSparseSpaceType::IsDistributedSpace());
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosGetRank, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    Epetra_MpiComm epetra_comm(raw_mpi_comm);
    KRATOS_EXPECT_EQ(r_comm.Rank(), TrilinosSparseSpaceType::GetRank(epetra_comm));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosIsNull, KratosTrilinosApplicationMPITestSuite)
{
    TrilinosSparseSpaceType::MatrixPointerType pA = nullptr;
    KRATOS_EXPECT_TRUE(TrilinosSparseSpaceType::IsNull(pA));
    // CreateEmptyMatrixPointer() intentionally returns nullptr (no-op empty pointer)
    pA = TrilinosSparseSpaceType::CreateEmptyMatrixPointer();
    KRATOS_EXPECT_TRUE(TrilinosSparseSpaceType::IsNull(pA));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosCreateEmptyPointer, KratosTrilinosApplicationMPITestSuite)
{
    auto p_map = TrilinosSparseSpaceType::CreateEmptyMapPointer();
    KRATOS_EXPECT_TRUE(TrilinosSparseSpaceType::IsNull(p_map));
    // CreateEmptyMatrixPointer/VectorPointer() return nullptr — IsNull = true
    auto p_A = TrilinosSparseSpaceType::CreateEmptyMatrixPointer();
    KRATOS_EXPECT_TRUE(TrilinosSparseSpaceType::IsNull(p_A));
    auto p_v = TrilinosSparseSpaceType::CreateEmptyVectorPointer();
    KRATOS_EXPECT_TRUE(TrilinosSparseSpaceType::IsNull(p_v));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosGetMap, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    auto vector = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size);
    auto& r_map = TrilinosSparseSpaceType::GetMap(vector);
    KRATOS_EXPECT_EQ(size, r_map.NumGlobalElements());
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosGetCommunicator, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    auto vector = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size);
    auto& r_trilinos_comm = TrilinosSparseSpaceType::GetCommunicator(vector);
    KRATOS_EXPECT_EQ(r_comm.Rank(), r_trilinos_comm.MyPID());
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosGlobalAssemble, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    auto vector = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size);
    auto matrix = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size);
    TrilinosSparseSpaceType::GlobalAssemble(vector);
    TrilinosSparseSpaceType::GlobalAssemble(matrix);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosBuildSystemStructure, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    Epetra_MpiComm epetra_comm(raw_mpi_comm);

    const int local_size = 2;
    const int first_my_id = r_comm.Rank() * local_size;
    const int system_size = 2 * r_comm.Size();

    std::vector<std::vector<int>> all_equation_ids = {{first_my_id, first_my_id + 1}};

    TrilinosSparseSpaceType::MatrixPointerType pA;
    TrilinosSparseSpaceType::VectorPointerType pb;
    TrilinosSparseSpaceType::VectorPointerType pDx;
    TrilinosSparseSpaceType::VectorPointerType pReactions;

    TrilinosSparseSpaceType::MapPointerType pMap;
    std::vector<int> local_ids(local_size);
    for (int i = 0; i < local_size; i++) local_ids[i] = first_my_id + i;
    pMap = Kratos::make_shared<Epetra_Map>(-1, local_size, local_ids.data(), 0, epetra_comm);

    TrilinosSparseSpaceType::BuildSystemStructure(
        epetra_comm, local_size, first_my_id, 5, all_equation_ids, all_equation_ids,
        pA, pb, pDx, pReactions, system_size, pMap
    );

    KRATOS_EXPECT_FALSE(TrilinosSparseSpaceType::IsNull(pA));
    KRATOS_EXPECT_FALSE(TrilinosSparseSpaceType::IsNull(pb));
    KRATOS_EXPECT_FALSE(TrilinosSparseSpaceType::IsNull(pDx));
    KRATOS_EXPECT_FALSE(TrilinosSparseSpaceType::IsNull(pReactions));
    KRATOS_EXPECT_EQ(static_cast<std::size_t>(system_size), TrilinosSparseSpaceType::Size1(*pA));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosBuildSystemStructureRowColumnBlocks, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    Epetra_MpiComm epetra_comm(raw_mpi_comm);

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
    std::vector<int> local_ids(local_size);
    for (int i = 0; i < local_size; i++) local_ids[i] = first_my_id + i;
    pMap = Kratos::make_shared<Epetra_Map>(-1, local_size, local_ids.data(), 0, epetra_comm);

    TrilinosSparseSpaceType::BuildSystemStructure(
        epetra_comm, local_size, first_my_id, 5, all_row_equation_ids, all_col_equation_ids,
        pA, pb, pDx, pReactions, system_size, pMap
    );

    KRATOS_EXPECT_FALSE(TrilinosSparseSpaceType::IsNull(pA));
    KRATOS_EXPECT_FALSE(TrilinosSparseSpaceType::IsNull(pb));
    KRATOS_EXPECT_FALSE(TrilinosSparseSpaceType::IsNull(pDx));
    KRATOS_EXPECT_FALSE(TrilinosSparseSpaceType::IsNull(pReactions));
    KRATOS_EXPECT_EQ(static_cast<std::size_t>(system_size), TrilinosSparseSpaceType::Size1(*pA));
    KRATOS_EXPECT_EQ(2 * r_comm.Size(), pA->NumGlobalNonzeros());
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosBuildConstraintsStructure, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    Epetra_MpiComm epetra_comm(raw_mpi_comm);

    const int local_size = 2;
    const int first_my_id = r_comm.Rank() * local_size;

    std::vector<std::vector<int>> all_slave_ids = {{first_my_id}};
    std::vector<std::vector<int>> all_master_ids = {{first_my_id + 1}};

    TrilinosSparseSpaceType::MatrixPointerType pT;
    TrilinosSparseSpaceType::VectorPointerType pConstantVector;

    TrilinosSparseSpaceType::MapPointerType pMap;
    std::vector<int> local_ids(local_size);
    for (int i = 0; i < local_size; i++) local_ids[i] = first_my_id + i;
    pMap = Kratos::make_shared<Epetra_Map>(-1, local_size, local_ids.data(), 0, epetra_comm);

    TrilinosSparseSpaceType::BuildConstraintsStructure(
        epetra_comm, local_size, first_my_id, 5, all_slave_ids, all_master_ids,
        pT, pConstantVector, pMap
    );

    KRATOS_EXPECT_FALSE(TrilinosSparseSpaceType::IsNull(pT));
    KRATOS_EXPECT_FALSE(TrilinosSparseSpaceType::IsNull(pConstantVector));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosManualFinalize, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    auto matrix = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size);
    // No-op for Epetra; ensure it doesn't throw
    TrilinosSparseSpaceType::ManualFinalize(matrix);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosGetCommunicatorMatrix, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    auto matrix = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size);
    auto& r_trilinos_comm = TrilinosSparseSpaceType::GetCommunicator(matrix);
    KRATOS_EXPECT_EQ(r_comm.Rank(), r_trilinos_comm.MyPID());
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosCreateVectorCopy, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    auto vector = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size);
    auto p_copy = TrilinosSparseSpaceType::CreateVectorCopy(vector);
    auto local_vector = TrilinosCPPTestUtilities::GenerateDummyLocalVector(size);
    TrilinosCPPTestUtilities::CheckSparseVectorFromLocalVector(*p_copy, local_vector);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosCreateEmptyPointerWithComm, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    TrilinosSparseSpaceType::CommunicatorType epetra_comm(raw_mpi_comm);

    auto pA = TrilinosSparseSpaceType::CreateEmptyMatrixPointer(epetra_comm);
    KRATOS_EXPECT_FALSE(TrilinosSparseSpaceType::IsNull(pA));
    KRATOS_EXPECT_EQ(0, static_cast<int>(TrilinosSparseSpaceType::Size1(*pA)));

    auto pv = TrilinosSparseSpaceType::CreateEmptyVectorPointer(epetra_comm);
    KRATOS_EXPECT_FALSE(TrilinosSparseSpaceType::IsNull(pv));
    KRATOS_EXPECT_EQ(0, static_cast<int>(TrilinosSparseSpaceType::Size(*pv)));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosCopyVector, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    auto vector_x = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size);
    // vector_y starts with offset=1 values, should be overwritten by Copy
    auto vector_y = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size, 1.0);
    TrilinosSparseSpaceType::Copy(vector_x, vector_y);
    auto local_vector = TrilinosCPPTestUtilities::GenerateDummyLocalVector(size);
    TrilinosCPPTestUtilities::CheckSparseVectorFromLocalVector(vector_y, local_vector);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosCopyMatrix, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    auto matrix_x = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size);
    // matrix_y starts with offset=1 values, should be overwritten by Copy
    auto matrix_y = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 1.0);
    TrilinosSparseSpaceType::Copy(matrix_x, matrix_y);
    auto local_matrix = TrilinosCPPTestUtilities::GenerateDummyLocalMatrix(size);
    TrilinosCPPTestUtilities::CheckSparseMatrixFromLocalMatrix(matrix_y, local_matrix);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosSetValue, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    TrilinosSparseSpaceType::CommunicatorType epetra_comm(raw_mpi_comm);
    TrilinosSparseSpaceType::MapType map(size, 0, epetra_comm);
    TrilinosVectorType vector(map);

    // Each rank sets its first owned global index to 42.0
    const int global_idx = r_comm.Rank() * 2;
    TrilinosSparseSpaceType::SetValue(vector, global_idx, 42.0);

    if (map.MyGID(global_idx)) {
        KRATOS_EXPECT_DOUBLE_EQ(42.0, vector[0][map.LID(global_idx)]);
    }
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosSetGlobalVecMethods, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    TrilinosSparseSpaceType::CommunicatorType epetra_comm(raw_mpi_comm);
    TrilinosSparseSpaceType::MapType map(size, 0, epetra_comm);
    TrilinosVectorType vector(map);
    vector.PutScalar(0.0);

    const int first_my_gid = r_comm.Rank() * 2;
    const int second_my_gid = first_my_gid + 1;

    TrilinosSparseSpaceType::SetGlobalVec(vector, first_my_gid, 42.0);
    if (map.MyGID(first_my_gid)) {
        KRATOS_EXPECT_DOUBLE_EQ(42.0, vector[0][map.LID(first_my_gid)]);
    }

    TrilinosSparseSpaceType::SetGlobalVecNoAssemble(vector, second_my_gid, 21.0);
    vector.GlobalAssemble(Insert, true);
    if (map.MyGID(second_my_gid)) {
        KRATOS_EXPECT_DOUBLE_EQ(21.0, vector[0][map.LID(second_my_gid)]);
    }
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosSetLocalVecMethods, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    TrilinosSparseSpaceType::CommunicatorType epetra_comm(raw_mpi_comm);
    TrilinosSparseSpaceType::MapType map(size, 0, epetra_comm);
    TrilinosVectorType vector(map);
    vector.PutScalar(0.0);

    TrilinosSparseSpaceType::SetLocalVec(vector, 0, 11.0);
    KRATOS_EXPECT_DOUBLE_EQ(11.0, vector[0][0]);

    TrilinosSparseSpaceType::SetLocalVecNoAssemble(vector, 1, 7.0);
    vector.GlobalAssemble(Insert, true);
    KRATOS_EXPECT_DOUBLE_EQ(7.0, vector[0][1]);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosSetGlobalMatMethods, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    const int first_my_gid = r_comm.Rank() * 2;
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    TrilinosSparseSpaceType::CommunicatorType epetra_comm(raw_mpi_comm);
    const TrilinosSparseSpaceType::MapType map(size, 0, epetra_comm);

    Epetra_FECrsGraph graph(Copy, map, 2);
    int col_gids[2] = {first_my_gid, first_my_gid + 1};
    for (int gid = first_my_gid; gid < first_my_gid + 2; ++gid) {
        graph.InsertGlobalIndices(1, &gid, 2, col_gids);
    }
    graph.GlobalAssemble();
    graph.FillComplete();
    TrilinosSparseMatrixType A(Copy, graph);
    A.PutScalar(0.0);

    TrilinosSparseSpaceType::SetGlobalMat(A, first_my_gid, first_my_gid, 5.0);

    int num_entries = 0;
    double* vals = nullptr;
    int* cols = nullptr;
    A.ExtractMyRowView(map.LID(first_my_gid), num_entries, vals, cols);
    for (int k = 0; k < num_entries; ++k) {
        if (A.ColMap().GID(cols[k]) == first_my_gid) {
            KRATOS_EXPECT_DOUBLE_EQ(5.0, vals[k]);
        }
    }

    TrilinosSparseSpaceType::SetGlobalMatNoAssemble(A, first_my_gid, first_my_gid + 1, 3.0);
    A.GlobalAssemble();

    A.ExtractMyRowView(map.LID(first_my_gid), num_entries, vals, cols);
    for (int k = 0; k < num_entries; ++k) {
        if (A.ColMap().GID(cols[k]) == first_my_gid + 1) {
            KRATOS_EXPECT_DOUBLE_EQ(3.0, vals[k]);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosSetLocalMatMethods, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    const int first_my_gid = r_comm.Rank() * 2;
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    TrilinosSparseSpaceType::CommunicatorType epetra_comm(raw_mpi_comm);
    const TrilinosSparseSpaceType::MapType map(size, 0, epetra_comm);

    Epetra_FECrsGraph graph(Copy, map, 2);
    int col_gids[2] = {first_my_gid, first_my_gid + 1};
    for (int gid = first_my_gid; gid < first_my_gid + 2; ++gid) {
        graph.InsertGlobalIndices(1, &gid, 2, col_gids);
    }
    graph.GlobalAssemble();
    graph.FillComplete();
    TrilinosSparseMatrixType A(Copy, graph);
    A.PutScalar(0.0);

    TrilinosSparseSpaceType::SetLocalMat(A, 0, 0, 9.0);

    int num_entries = 0;
    double* vals = nullptr;
    int* cols = nullptr;
    A.ExtractMyRowView(0, num_entries, vals, cols);
    for (int k = 0; k < num_entries; ++k) {
        if (cols[k] == 0) {
            KRATOS_EXPECT_DOUBLE_EQ(9.0, vals[k]);
        }
    }

    TrilinosSparseSpaceType::SetLocalMatNoAssemble(A, 0, 1, 4.0);
    A.GlobalAssemble();

    A.ExtractMyRowView(0, num_entries, vals, cols);
    for (int k = 0; k < num_entries; ++k) {
        if (cols[k] == 1) {
            KRATOS_EXPECT_DOUBLE_EQ(4.0, vals[k]);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosGetValue, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    auto vector = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size);
    // Each rank only owns its 2 local indices
    const int first_my_id = r_comm.Rank() * 2;
    KRATOS_EXPECT_DOUBLE_EQ(static_cast<double>(first_my_id),
        TrilinosSparseSpaceType::GetValue(vector, first_my_id));
    KRATOS_EXPECT_DOUBLE_EQ(static_cast<double>(first_my_id + 1),
        TrilinosSparseSpaceType::GetValue(vector, first_my_id + 1));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosGatherValues, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    auto vector = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size);
    // Gather the first two global values on every process
    std::vector<int> indices = {0, 1};
    std::vector<double> values(2, 0.0);
    TrilinosSparseSpaceType::GatherValues(vector, indices, values.data());
    KRATOS_EXPECT_DOUBLE_EQ(0.0, values[0]);
    KRATOS_EXPECT_DOUBLE_EQ(1.0, values[1]);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosClear, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    TrilinosSparseSpaceType::CommunicatorType epetra_comm(raw_mpi_comm);

    auto pA = TrilinosSparseSpaceType::CreateEmptyMatrixPointer(epetra_comm);
    KRATOS_EXPECT_FALSE(TrilinosSparseSpaceType::IsNull(pA));
    TrilinosSparseSpaceType::Clear(pA);
    KRATOS_EXPECT_FALSE(TrilinosSparseSpaceType::IsNull(pA)); // still valid, just empty
    KRATOS_EXPECT_EQ(0, static_cast<int>(TrilinosSparseSpaceType::Size1(*pA)));

    auto pv = TrilinosSparseSpaceType::CreateEmptyVectorPointer(epetra_comm);
    KRATOS_EXPECT_FALSE(TrilinosSparseSpaceType::IsNull(pv));
    TrilinosSparseSpaceType::Clear(pv);
    KRATOS_EXPECT_FALSE(TrilinosSparseSpaceType::IsNull(pv)); // still valid, just empty
    KRATOS_EXPECT_EQ(0, static_cast<int>(TrilinosSparseSpaceType::Size(*pv)));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosResizeVectorPointer, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    TrilinosSparseSpaceType::CommunicatorType epetra_comm(raw_mpi_comm);

    auto pv = TrilinosSparseSpaceType::CreateEmptyVectorPointer(epetra_comm);
    KRATOS_EXPECT_EQ(0, static_cast<int>(TrilinosSparseSpaceType::Size(*pv)));
    TrilinosSparseSpaceType::Resize(pv, size);
    KRATOS_EXPECT_EQ(static_cast<std::size_t>(size), TrilinosSparseSpaceType::Size(*pv));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosFastestDirectSolverList, KratosTrilinosApplicationMPITestSuite)
{
    const auto solvers = TrilinosSparseSpaceType::FastestDirectSolverList();
    KRATOS_EXPECT_GT(solvers.size(), std::size_t(0));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosGetDiagonalNorm, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 12;
    auto matrix = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 1.0);
    // Diagonal values are 1, 2, ..., 12 → norm = sqrt(1² + 2² + ... + 12²)
    double ref = 0.0;
    for (int i = 1; i <= size; ++i) ref += static_cast<double>(i) * static_cast<double>(i);
    ref = std::sqrt(ref);
    KRATOS_EXPECT_NEAR(TrilinosSparseSpaceType::GetDiagonalNorm(matrix), ref, 1.0e-6);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosGetMaxDiagonal, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 12;
    auto matrix = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(12.0, TrilinosSparseSpaceType::GetMaxDiagonal(matrix));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosAssembleLHS, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    const int first_my_id = r_comm.Rank() * 2;
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    TrilinosSparseSpaceType::CommunicatorType epetra_comm(raw_mpi_comm);
    const TrilinosSparseSpaceType::MapType map(size, 0, epetra_comm);

    // Build FE graph with full 2x2 dense block per rank
    Epetra_FECrsGraph graph(Copy, map, 2);
    int col_gids[2] = {first_my_id, first_my_id + 1};
    for (int gid = first_my_id; gid < first_my_id + 2; ++gid) {
        graph.InsertGlobalIndices(1, &gid, 2, col_gids);
    }
    graph.GlobalAssemble();
    graph.FillComplete();
    TrilinosSparseMatrixType A(Copy, graph);

    // Assemble 2x2 local diagonal-like matrix
    Matrix local_mat(2, 2, 0.0);
    local_mat(0, 0) = 5.0;
    local_mat(1, 1) = 3.0;
    std::vector<std::size_t> eq_ids = {static_cast<std::size_t>(first_my_id),
                                        static_cast<std::size_t>(first_my_id + 1)};
    TrilinosSparseSpaceType::AssembleLHS(A, local_mat, eq_ids);
    A.GlobalAssemble();

    if (map.MyGID(first_my_id)) {
        int num_entries;
        double* vals;
        int* cols;
        A.ExtractMyRowView(map.LID(first_my_id), num_entries, vals, cols);
        for (int k = 0; k < num_entries; ++k) {
            if (A.ColMap().GID(cols[k]) == first_my_id) {
                KRATOS_EXPECT_DOUBLE_EQ(5.0, vals[k]);
            }
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosAssembleRHS, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    const int first_my_id = r_comm.Rank() * 2;
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    TrilinosSparseSpaceType::CommunicatorType epetra_comm(raw_mpi_comm);
    TrilinosSparseSpaceType::MapType map(size, 0, epetra_comm);
    TrilinosVectorType b(map);

    Vector local_vec(2);
    local_vec[0] = 7.0;
    local_vec[1] = 3.0;
    std::vector<std::size_t> eq_ids = {static_cast<std::size_t>(first_my_id),
                                        static_cast<std::size_t>(first_my_id + 1)};
    TrilinosSparseSpaceType::AssembleRHS(b, local_vec, eq_ids);
    b.GlobalAssemble();

    if (map.MyGID(first_my_id)) {
        KRATOS_EXPECT_DOUBLE_EQ(7.0, b[0][map.LID(first_my_id)]);
    }
    if (map.MyGID(first_my_id + 1)) {
        KRATOS_EXPECT_DOUBLE_EQ(3.0, b[0][map.LID(first_my_id + 1)]);
    }
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosCreateDofUpdater, KratosTrilinosApplicationMPITestSuite)
{
    auto dof_updater = TrilinosSparseSpaceType::CreateDofUpdater();
    KRATOS_EXPECT_NE(nullptr, dof_updater.get());
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosGetColumnThrows, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    TrilinosSparseSpaceType::CommunicatorType epetra_comm(raw_mpi_comm);
    TrilinosSparseSpaceType::MapType map(size, 0, epetra_comm);
    auto matrix = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size);
    TrilinosVectorType col(map);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        TrilinosSparseSpaceType::GetColumn(0, matrix, col),
        "GetColumn method is not currently implemented"
    );
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosResizeThrows, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 2 * r_comm.Size();
    auto matrix = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size);
    auto vector = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        TrilinosSparseSpaceType::Resize(matrix, size, size),
        "Resize is not defined for Trilinos Sparse Matrix"
    );
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        TrilinosSparseSpaceType::Resize(vector, size),
        "Resize is not defined for a reference to Trilinos Vector"
    );
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosCreateEmptyMapPointer, KratosTrilinosApplicationMPITestSuite)
{
    // CreateEmptyMapPointer always returns nullptr
    auto p_map = TrilinosSparseSpaceType::CreateEmptyMapPointer();
    KRATOS_EXPECT_TRUE(p_map == nullptr);
    KRATOS_EXPECT_TRUE(TrilinosSparseSpaceType::IsNull(p_map));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosGetMinDiagonal, KratosTrilinosApplicationMPITestSuite)
{
    // Diagonal of the dummy matrix (scale=1) is i+1 for global row i, so min = 1.0
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 12;
    auto matrix = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(1.0, TrilinosSparseSpaceType::GetMinDiagonal(matrix));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosGetAveragevalueDiagonal, KratosTrilinosApplicationMPITestSuite)
{
    // Diagonal entries are 1, 2, ..., 12 → average = (1+12)/2 = 6.5
    const auto& r_comm = Testing::GetDefaultDataCommunicator();
    const int size = 12;
    auto matrix = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 1.0);
    const double expected_avg = 0.5 * (1.0 + 12.0);
    KRATOS_EXPECT_DOUBLE_EQ(expected_avg, TrilinosSparseSpaceType::GetAveragevalueDiagonal(matrix));
}

} // namespace Kratos::Testing