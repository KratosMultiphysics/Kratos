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
#include "tests/cpp_tests/trilinos_fast_suite.h"
#include "trilinos_space_experimental.h"
#include "containers/model.h"
#include "mpi/includes/mpi_data_communicator.h"
#include "trilinos_cpp_test_experimental_utilities.h"
#include "utilities/math_utils.h"

namespace Kratos::Testing
{

/// Basic definitions
using TrilinosSparseMatrix = Tpetra::FECrsMatrix<>;
using TrilinosSparseVector = Tpetra::Vector<>;
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
    auto matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 0.0, true);
    auto local_matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalMatrix(size, 0.0, true);

    // Check
    KRATOS_EXPECT_DOUBLE_EQ(TrilinosLocalSpaceType::TwoNorm(local_matrix), TrilinosSparseSpaceType::TwoNorm(*matrix));
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalMultMatrixVector, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 2 * r_comm.Size();
    auto matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 0.0, true);
    auto local_matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalMatrix(size, 0.0, true);
    auto vector = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size, 0.0);
    auto local_vector = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalVector(size, 0.0);

    // Tpetra communicator
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    TrilinosSparseSpaceType::CommunicatorPointerType tpetra_comm = Teuchos::rcp(new TrilinosSparseSpaceType::CommunicatorType(raw_mpi_comm));

    // Create a map
    const int global_elems = 0;
    TrilinosSparseSpaceType::MapPointerType map = Teuchos::rcp(new TrilinosSparseSpaceType::MapType(global_elems, 0, tpetra_comm));

    // Create an Tpetra_Vector
    TrilinosVectorType mult(map);

    // Solution
    TrilinosSparseSpaceType::Mult(*matrix, *vector, mult);

    // Check
    const TrilinosLocalVectorType multiply_reference = prod(local_matrix, local_vector);
    TrilinosCPPTestExperimentalUtilities::CheckSparseVectorFromLocalVector(mult, multiply_reference);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalMultMatrixMatrix, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 2 * r_comm.Size();
    auto matrix_1 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 0.0, true);
    auto local_matrix_1 = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalMatrix(size, 0.0, true);
    auto matrix_2 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 1.0, true);
    auto local_matrix_2 = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalMatrix(size, 1.0, true);

    // Tpetra communicator
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    TrilinosSparseSpaceType::CommunicatorPointerType tpetra_comm = Teuchos::rcp(new TrilinosSparseSpaceType::CommunicatorType(raw_mpi_comm));

    // Create a map
    const int global_elems = 0;
    TrilinosSparseSpaceType::MapPointerType map = Teuchos::rcp(new TrilinosSparseSpaceType::MapType(global_elems, 0, tpetra_comm));

    // Create an Tpetra_Matrix
    TrilinosSparseSpaceType::GraphPointerType graph = Teuchos::rcp(new TrilinosSparseSpaceType::GraphType(map, map, 0));
    TrilinosSparseMatrixType mult(graph);

    // Solution
    TrilinosSparseSpaceType::Mult(*matrix_1, *matrix_2, mult);

    // Check
    const TrilinosLocalMatrixType multiply_reference = prod(local_matrix_1, local_matrix_2);
    TrilinosCPPTestExperimentalUtilities::CheckSparseMatrixFromLocalMatrix(mult, multiply_reference);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalTransposeMultMatrixVector, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 2 * r_comm.Size();
    auto matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 0.0, true);
    auto local_matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalMatrix(size, 0.0, true);
    auto vector = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size, 0.0);
    auto local_vector = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalVector(size, 0.0);

    // Tpetra communicator
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    TrilinosSparseSpaceType::CommunicatorPointerType tpetra_comm = Teuchos::rcp(new TrilinosSparseSpaceType::CommunicatorType(raw_mpi_comm));

    // Create a map
    const int global_elems = 0;
    TrilinosSparseSpaceType::MapPointerType map = Teuchos::rcp(new TrilinosSparseSpaceType::MapType(global_elems, 0, tpetra_comm));

    // Create an Tpetra_Vector
    TrilinosVectorType mult(map);

    // Solution
    TrilinosSparseSpaceType::TransposeMult(*matrix, *vector, mult);

    // Check
    const TrilinosLocalVectorType multiply_reference = prod(trans(local_matrix), local_vector);
    TrilinosCPPTestExperimentalUtilities::CheckSparseVectorFromLocalVector(mult, multiply_reference);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalTransposeMultMatrixMatrix, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 2 * r_comm.Size();
    auto matrix_1 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 0.0, true);
    auto local_matrix_1 = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalMatrix(size, 0.0, true);
    auto matrix_2 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 1.0, true);
    auto local_matrix_2 = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalMatrix(size, 1.0, true);

    // Tpetra communicator
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    TrilinosSparseSpaceType::CommunicatorPointerType tpetra_comm = Teuchos::rcp(new TrilinosSparseSpaceType::CommunicatorType(raw_mpi_comm));

    // Create a map
    const int global_elems = 0;
    TrilinosSparseSpaceType::MapPointerType map = Teuchos::rcp(new TrilinosSparseSpaceType::MapType(global_elems, 0, tpetra_comm));

    // Create an Tpetra_Matrix
    TrilinosSparseSpaceType::GraphPointerType graph = Teuchos::rcp(new TrilinosSparseSpaceType::GraphType(map, map, 0));
    TrilinosSparseMatrixType mult(graph);

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
    auto matrix_1 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 0.0, true);
    auto local_matrix_1 = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalMatrix(size, 0.0, true);
    auto matrix_2 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 1.0, true);
    auto local_matrix_2 = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalMatrix(size, 1.0, true);

    // Tpetra communicator
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    TrilinosSparseSpaceType::CommunicatorPointerType tpetra_comm = Teuchos::rcp(new TrilinosSparseSpaceType::CommunicatorType(raw_mpi_comm));

    // Create a map
    const int global_elems = 0;
    TrilinosSparseSpaceType::MapPointerType map = Teuchos::rcp(new TrilinosSparseSpaceType::MapType(global_elems, 0, tpetra_comm));

    // Create an Tpetra_Matrix
    TrilinosSparseSpaceType::GraphPointerType graph = Teuchos::rcp(new TrilinosSparseSpaceType::GraphType(map, map, 0));
    TrilinosSparseMatrixType mult(graph);

    // Solution
    TrilinosSparseSpaceType::BtDBProductOperation(mult, *matrix_1, *matrix_2);

    // Check
    TrilinosLocalMatrixType multiply_reference;
    MathUtils<double>::BtDBProductOperation(multiply_reference, local_matrix_1, local_matrix_2);
    TrilinosCPPTestExperimentalUtilities::CheckSparseMatrixFromLocalMatrix(mult, multiply_reference);

    // Non zero matrix
    auto second_mult = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 0.0, true);

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
    auto matrix_1 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 0.0, true);
    auto local_matrix_1 = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalMatrix(size, 0.0, true);
    auto matrix_2 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 1.0, true);
    auto local_matrix_2 = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalMatrix(size, 1.0, true);

    // Tpetra communicator
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    TrilinosSparseSpaceType::CommunicatorPointerType tpetra_comm = Teuchos::rcp(new TrilinosSparseSpaceType::CommunicatorType(raw_mpi_comm));

    // Create a map
    const int global_elems = 0;
    TrilinosSparseSpaceType::MapPointerType map = Teuchos::rcp(new TrilinosSparseSpaceType::MapType(global_elems, 0, tpetra_comm));

    // Create an Tpetra_Matrix
    TrilinosSparseSpaceType::GraphPointerType graph = Teuchos::rcp(new TrilinosSparseSpaceType::GraphType(map, map, 0));
    TrilinosSparseMatrixType mult(graph);

    // Solution
    TrilinosSparseSpaceType::BDBtProductOperation(mult, *matrix_1, *matrix_2);

    // Check
    TrilinosLocalMatrixType multiply_reference;
    MathUtils<double>::BDBtProductOperation(multiply_reference, local_matrix_1, local_matrix_2);
    TrilinosCPPTestExperimentalUtilities::CheckSparseMatrixFromLocalMatrix(mult, multiply_reference);

    // Non zero matrix
    auto second_mult = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 0.0, true);

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
    auto matrix_1 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 0.0, true);
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
    auto matrix_2 = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size, 0.0, true);
    auto local_matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalMatrix(size, 0.0, true);

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
    TrilinosVectorType vector12(map);

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

} // namespace Kratos::Testing