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
#include "testing/testing.h"
#include "trilinos_space.h"
#include "containers/model.h"
#include "mpi/includes/mpi_data_communicator.h"
#include "trilinos_cpp_test_utilities.h"
#include "utilities/math_utils.h"

namespace Kratos::Testing
{

/// Basic definitions
using TrilinosSparseSpaceType = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>;
using TrilinosLocalSpaceType = UblasSpace<double, Matrix, Vector>;

using TrilinosSparseMatrixType = TrilinosSparseSpaceType::MatrixType;
using TrilinosVectorType = TrilinosSparseSpaceType::VectorType;

using TrilinosLocalMatrixType = TrilinosLocalSpaceType::MatrixType;
using TrilinosLocalVectorType = TrilinosLocalSpaceType::VectorType;

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosSizeVector, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto vector = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size);

    // Check
    KRATOS_EXPECT_EQ(static_cast<std::size_t>(size), TrilinosSparseSpaceType::Size(vector));
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosSizeMatrix, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto vector = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size);

    // Check
    KRATOS_EXPECT_EQ(static_cast<std::size_t>(size), TrilinosSparseSpaceType::Size1(vector));
    KRATOS_EXPECT_EQ(static_cast<std::size_t>(size), TrilinosSparseSpaceType::Size2(vector));
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosDotProduct, KratosTrilinosApplicationMPITestSuite)
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

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosMaxMin, KratosTrilinosApplicationMPITestSuite)
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

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosTwoNormVector, KratosTrilinosApplicationMPITestSuite)
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

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosTwoNormMatrix1, KratosTrilinosApplicationMPITestSuite)
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

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosTwoNormMatrix2, KratosTrilinosApplicationMPITestSuite)
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

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosMultMatrixVector, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 2 * r_comm.Size();
    auto matrix = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 0.0, true);
    auto local_matrix = TrilinosCPPTestUtilities::GenerateDummyLocalMatrix(size, 0.0, true);
    auto vector = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size, 0.0);
    auto local_vector = TrilinosCPPTestUtilities::GenerateDummyLocalVector(size, 0.0);

    // Epetra coomunicator
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    Epetra_MpiComm epetra_comm(raw_mpi_comm);

    // Create a map
    Epetra_Map Map(size,0,epetra_comm);

    // Create an Epetra_Vector
    TrilinosVectorType mult(Map);

    // Solution
    TrilinosSparseSpaceType::Mult(matrix, vector, mult);

    // Check
    const TrilinosLocalVectorType multiply_reference = prod(local_matrix, local_vector);
    TrilinosCPPTestUtilities::CheckSparseVectorFromLocalVector(mult, multiply_reference);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosMultMatrixMatrix, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 2 * r_comm.Size();
    auto matrix_1 = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 0.0, true);
    auto local_matrix_1 = TrilinosCPPTestUtilities::GenerateDummyLocalMatrix(size, 0.0, true);
    auto matrix_2 = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 1.0, true);
    auto local_matrix_2 = TrilinosCPPTestUtilities::GenerateDummyLocalMatrix(size, 1.0, true);

    // Epetra coomunicator
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    Epetra_MpiComm epetra_comm(raw_mpi_comm);

    // Create a map
    Epetra_Map Map(size,0,epetra_comm);

    // Create an Epetra_Matrix
    std::vector<int> NumNz;
    TrilinosSparseMatrixType mult(Copy, Map, NumNz.data());

    // Solution
    TrilinosSparseSpaceType::Mult(matrix_1, matrix_2, mult);

    // Check
    const TrilinosLocalMatrixType multiply_reference = prod(local_matrix_1, local_matrix_2);
    TrilinosCPPTestUtilities::CheckSparseMatrixFromLocalMatrix(mult, multiply_reference);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosTransposeMultMatrixVector, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 2 * r_comm.Size();
    auto matrix = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 0.0, true);
    auto local_matrix = TrilinosCPPTestUtilities::GenerateDummyLocalMatrix(size, 0.0, true);
    auto vector = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size, 0.0);
    auto local_vector = TrilinosCPPTestUtilities::GenerateDummyLocalVector(size, 0.0);

    // Epetra coomunicator
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    Epetra_MpiComm epetra_comm(raw_mpi_comm);

    // Create a map
    Epetra_Map Map(size,0,epetra_comm);

    // Create an Epetra_Vector
    TrilinosVectorType mult(Map);

    // Solution
    TrilinosSparseSpaceType::TransposeMult(matrix, vector, mult);

    // Check
    const TrilinosLocalVectorType multiply_reference = prod(trans(local_matrix), local_vector);
    TrilinosCPPTestUtilities::CheckSparseVectorFromLocalVector(mult, multiply_reference);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosTransposeMultMatrixMatrix, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 2 * r_comm.Size();
    auto matrix_1 = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 0.0, true);
    auto local_matrix_1 = TrilinosCPPTestUtilities::GenerateDummyLocalMatrix(size, 0.0, true);
    auto matrix_2 = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 1.0, true);
    auto local_matrix_2 = TrilinosCPPTestUtilities::GenerateDummyLocalMatrix(size, 1.0, true);

    // Epetra coomunicator
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    Epetra_MpiComm epetra_comm(raw_mpi_comm);

    // Create a map
    Epetra_Map Map(size,0,epetra_comm);

    // Create an Epetra_Matrix
    std::vector<int> NumNz;
    TrilinosSparseMatrixType mult(Copy, Map, NumNz.data());

    // Solution
    TrilinosSparseSpaceType::TransposeMult(matrix_1, matrix_2, mult, {true, false});

    // Check
    const TrilinosLocalMatrixType multiply_reference = prod(trans(local_matrix_1), local_matrix_2);
    TrilinosCPPTestUtilities::CheckSparseMatrixFromLocalMatrix(mult, multiply_reference);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosBtDBProductOperation, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 2 * r_comm.Size();
    auto matrix_1 = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 0.0, true);
    auto local_matrix_1 = TrilinosCPPTestUtilities::GenerateDummyLocalMatrix(size, 0.0, true);
    auto matrix_2 = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 1.0, true);
    auto local_matrix_2 = TrilinosCPPTestUtilities::GenerateDummyLocalMatrix(size, 1.0, true);

    // Epetra coomunicator
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    Epetra_MpiComm epetra_comm(raw_mpi_comm);

    // Create a map
    Epetra_Map Map(size,0,epetra_comm);

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
// KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosBtDBProductOperationRealCase, KratosTrilinosApplicationMPITestSuite)
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
//     const Epetra_Map* p_map = &(A.RowMap());

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

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosBDBtProductOperation, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 2 * r_comm.Size();
    auto matrix_1 = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 0.0, true);
    auto local_matrix_1 = TrilinosCPPTestUtilities::GenerateDummyLocalMatrix(size, 0.0, true);
    auto matrix_2 = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 1.0, true);
    auto local_matrix_2 = TrilinosCPPTestUtilities::GenerateDummyLocalMatrix(size, 1.0, true);

    // Epetra coomunicator
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    Epetra_MpiComm epetra_comm(raw_mpi_comm);

    // Create a map
    Epetra_Map Map(size,0,epetra_comm);

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

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosInplaceMult, KratosTrilinosApplicationMPITestSuite)
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

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosAssign, KratosTrilinosApplicationMPITestSuite)
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

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosUnaliasedAdd, KratosTrilinosApplicationMPITestSuite)
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

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosScaleAndAdd1, KratosTrilinosApplicationMPITestSuite)
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

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosScaleAndAdd2, KratosTrilinosApplicationMPITestSuite)
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

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosSet, KratosTrilinosApplicationMPITestSuite)
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

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosSetToZeroMatrix, KratosTrilinosApplicationMPITestSuite)
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

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosSetToZeroVector, KratosTrilinosApplicationMPITestSuite)
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

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosCopyMatrixValues, KratosTrilinosApplicationMPITestSuite)
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

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosCombineMatricesGraphs, KratosTrilinosApplicationMPITestSuite)
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

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosCheckAndCorrectZeroDiagonalValues, KratosTrilinosApplicationMPITestSuite)
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
    Epetra_MpiComm epetra_comm(raw_mpi_comm);

    // Dummy vector
    Epetra_Map map(size, 0, epetra_comm);
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

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosIsDistributed, KratosTrilinosApplicationMPITestSuite)
{
    KRATOS_EXPECT_TRUE(TrilinosSparseSpaceType::IsDistributed());
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosGetScaleNorm, KratosTrilinosApplicationMPITestSuite)
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

} // namespace Kratos::Testing