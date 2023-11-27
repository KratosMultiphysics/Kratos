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
#include "mpi/includes/mpi_data_communicator.h"
#include "../trilinos_cpp_test_utilities.h"
#include "custom_utilities/trilinos_assembling_utilities.h"

namespace Kratos::Testing
{
/// Basic definitions
using TrilinosSparseSpaceType = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>;
using TrilinosLocalSpaceType = UblasSpace<double, Matrix, Vector>;

using TrilinosSparseMatrixType = TrilinosSparseSpaceType::MatrixType;
using TrilinosVectorType = TrilinosSparseSpaceType::VectorType;

using TrilinosLocalMatrixType = TrilinosLocalSpaceType::MatrixType;
using TrilinosLocalVectorType = TrilinosLocalSpaceType::VectorType;


KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosVectorSetValue, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    const int rank = r_comm.Rank();
    auto vector = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size);
    TrilinosLocalVectorType local_vector = ZeroVector(size);
    for (int i = 0; i < size; ++i) {
        local_vector[i] = 1.0;
    }

    // Solution global
    TrilinosSparseSpaceType::SetToZero(vector);
    for (int i = 0; i < 2; ++i) {
        TrilinosAssemblingUtilities::SetGlobalValue(vector, rank * 2 + i, 1.0);
    }

    // Check global
    TrilinosCPPTestUtilities::CheckSparseVectorFromLocalVector(vector, local_vector);

    // Solution local
    TrilinosSparseSpaceType::SetToZero(vector);
    for (int i = 0; i < 2; ++i) {
        TrilinosAssemblingUtilities::SetLocalValue(vector, i, 1.0);
    }

    // Check local
    TrilinosCPPTestUtilities::CheckSparseVectorFromLocalVector(vector, local_vector);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosVectorSetGlobalValueWithoutGlobalAssembly, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    const int rank = r_comm.Rank();
    auto vector = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size);
    TrilinosLocalVectorType local_vector = ZeroVector(size);
    for (int i = 0; i < size; ++i) {
        local_vector[i] = 1.0;
    }

    // Solution global
    TrilinosSparseSpaceType::SetToZero(vector);
    for (int i = 0; i < 2; ++i) {
        TrilinosAssemblingUtilities::SetGlobalValueWithoutGlobalAssembly(vector, rank * 2 + i, 1.0);
    }
    vector.GlobalAssemble();

    // Check global
    TrilinosCPPTestUtilities::CheckSparseVectorFromLocalVector(vector, local_vector);

    // Solution local
    TrilinosSparseSpaceType::SetToZero(vector);
    for (int i = 0; i < 2; ++i) {
        TrilinosAssemblingUtilities::SetLocalValueWithoutGlobalAssembly(vector, i, 1.0);
    }
    vector.GlobalAssemble();

    // Check local
    TrilinosCPPTestUtilities::CheckSparseVectorFromLocalVector(vector, local_vector);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosMatrixSetValue, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    const int rank = r_comm.Rank();
    auto matrix = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size);
    TrilinosLocalMatrixType local_matrix = ZeroMatrix(size, size);
    for (int i = 0; i < size; ++i) {
        local_matrix(i, i) = 1.0;
    }

    // Solution global
    for (int i = 0; i < 2; ++i) {
        TrilinosAssemblingUtilities::SetGlobalValue(matrix, rank * 2 + i, rank * 2 + i, 1.0);
    }

    // Check global
    TrilinosCPPTestUtilities::CheckSparseMatrixFromLocalMatrix(matrix, local_matrix);

    // Solution local
    TrilinosSparseSpaceType::SetToZero(matrix);
    for (int i = 0; i < 2; ++i) {
        TrilinosAssemblingUtilities::SetLocalValue(matrix, i, i, 1.0);
    }

    // Check local
    TrilinosCPPTestUtilities::CheckSparseMatrixFromLocalMatrix(matrix, local_matrix);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosMatrixSetGlobalValueWithoutGlobalAssembly, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    const int rank = r_comm.Rank();
    auto matrix = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size);
    TrilinosLocalMatrixType local_matrix = ZeroMatrix(size, size);
    for (int i = 0; i < size; ++i) {
        local_matrix(i, i) = 1.0;
    }

    // Solution global
    for (int i = 0; i < 2; ++i) {
        TrilinosAssemblingUtilities::SetGlobalValueWithoutGlobalAssembly(matrix, rank * 2 + i, rank * 2 + i, 1.0);
    }
    matrix.GlobalAssemble();

    // Check global
    TrilinosCPPTestUtilities::CheckSparseMatrixFromLocalMatrix(matrix, local_matrix);

    // Solution local
    TrilinosSparseSpaceType::SetToZero(matrix);
    for (int i = 0; i < 2; ++i) {
        TrilinosAssemblingUtilities::SetLocalValueWithoutGlobalAssembly(matrix, i, i, 1.0);
    }
    matrix.GlobalAssemble();

    // Check local
    TrilinosCPPTestUtilities::CheckSparseMatrixFromLocalMatrix(matrix, local_matrix);
}

} // namespace Kratos::Testing