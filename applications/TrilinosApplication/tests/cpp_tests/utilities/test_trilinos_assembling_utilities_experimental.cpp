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
#include "mpi/includes/mpi_data_communicator.h"
#include "../trilinos_cpp_test_experimental_utilities.h"
#include "custom_utilities/trilinos_assembling_utilities.h"
#include "trilinos_space_experimental.h"

namespace Kratos::Testing
{
/// Basic definitions
using TrilinosSparseSpaceType = TrilinosSpaceExperimental<Tpetra::FECrsMatrix<>, Tpetra::FEMultiVector<>>;
using TrilinosLocalSpaceType = UblasSpace<double, Matrix, Vector>;

using TrilinosSparseMatrixType = TrilinosSparseSpaceType::MatrixType;
using TrilinosVectorType = TrilinosSparseSpaceType::VectorType;

using TrilinosLocalMatrixType = TrilinosLocalSpaceType::MatrixType;
using TrilinosLocalVectorType = TrilinosLocalSpaceType::VectorType;


KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalVectorSetValue, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    const int rank = r_comm.Rank();
    auto vector = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size);
    TrilinosLocalVectorType local_vector = ZeroVector(size);
    for (int i = 0; i < size; ++i) {
        local_vector[i] = 1.0;
    }

    // Solution global
    TrilinosSparseSpaceType::SetToZero(*vector);
    for (int i = 0; i < 2; ++i) {
        TrilinosAssemblingUtilities<TrilinosSparseSpaceType>::SetGlobalValue(*vector, rank * 2 + i, 1.0);
    }

    // Check global
    TrilinosCPPTestExperimentalUtilities::CheckSparseVectorFromLocalVector(*vector, local_vector);

    // Solution local
    TrilinosSparseSpaceType::SetToZero(*vector);
    for (int i = 0; i < 2; ++i) {
        TrilinosAssemblingUtilities<TrilinosSparseSpaceType>::SetLocalValue(*vector, i, 1.0);
    }

    // Check local
    TrilinosCPPTestExperimentalUtilities::CheckSparseVectorFromLocalVector(*vector, local_vector);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalVectorSetGlobalValueWithoutGlobalAssembly, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    const int rank = r_comm.Rank();
    auto vector = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(r_comm, size);
    TrilinosLocalVectorType local_vector = ZeroVector(size);
    for (int i = 0; i < size; ++i) {
        local_vector[i] = 1.0;
    }

    // Solution global
    TrilinosSparseSpaceType::SetToZero(*vector);
    for (int i = 0; i < 2; ++i) {
        TrilinosAssemblingUtilities<TrilinosSparseSpaceType>::SetGlobalValueWithoutGlobalAssembly(*vector, rank * 2 + i, 1.0);
    }

    // Check global
    TrilinosCPPTestExperimentalUtilities::CheckSparseVectorFromLocalVector(*vector, local_vector);

    // Solution local
    TrilinosSparseSpaceType::SetToZero(*vector);
    for (int i = 0; i < 2; ++i) {
        TrilinosAssemblingUtilities<TrilinosSparseSpaceType>::SetLocalValueWithoutGlobalAssembly(*vector, i, 1.0);
    }

    // Check local
    TrilinosCPPTestExperimentalUtilities::CheckSparseVectorFromLocalVector(*vector, local_vector);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalMatrixSetValue, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    const int rank = r_comm.Rank();
    auto matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size);
    TrilinosLocalMatrixType local_matrix = ZeroMatrix(size, size);
    for (int i = 0; i < size; ++i) {
        local_matrix(i, i) = 1.0;
    }

    // Solution global
    if (!matrix->isFillActive()) matrix->resumeFill();
    matrix->beginAssembly();
    for (int i = 0; i < 2; ++i) {
        TrilinosAssemblingUtilities<TrilinosSparseSpaceType>::SetGlobalValue(*matrix, rank * 2 + i, rank * 2 + i, 1.0);
    }
    TrilinosSparseSpaceType::GlobalAssemble(*matrix);

    // Check global
    TrilinosCPPTestExperimentalUtilities::CheckSparseMatrixFromLocalMatrix(*matrix, local_matrix);

    // Solution local
    TrilinosSparseSpaceType::SetToZero(*matrix);
    if (!matrix->isFillActive()) matrix->resumeFill();
    matrix->beginAssembly();
    for (int i = 0; i < 2; ++i) {
        TrilinosAssemblingUtilities<TrilinosSparseSpaceType>::SetLocalValue(*matrix, i, i, 1.0);
    }
    TrilinosSparseSpaceType::GlobalAssemble(*matrix);

    // Check local
    TrilinosCPPTestExperimentalUtilities::CheckSparseMatrixFromLocalMatrix(*matrix, local_matrix);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalMatrixSetGlobalValueWithoutGlobalAssembly, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    const int rank = r_comm.Rank();
    auto matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size);
    TrilinosLocalMatrixType local_matrix = ZeroMatrix(size, size);
    for (int i = 0; i < size; ++i) {
        local_matrix(i, i) = 1.0;
    }

    // Solution global
    if (!matrix->isFillActive()) matrix->resumeFill();
    matrix->beginAssembly();
    for (int i = 0; i < 2; ++i) {
        TrilinosAssemblingUtilities<TrilinosSparseSpaceType>::SetGlobalValue(*matrix, rank * 2 + i, rank * 2 + i, 1.0);
    }
    TrilinosSparseSpaceType::GlobalAssemble(*matrix);

    // Check global
    TrilinosCPPTestExperimentalUtilities::CheckSparseMatrixFromLocalMatrix(*matrix, local_matrix);

    // Solution local
    TrilinosSparseSpaceType::SetToZero(*matrix);
    if (!matrix->isFillActive()) matrix->resumeFill();
    matrix->beginAssembly();
    for (int i = 0; i < 2; ++i) {
        TrilinosAssemblingUtilities<TrilinosSparseSpaceType>::SetLocalValue(*matrix, i, i, 1.0);
    }
    TrilinosSparseSpaceType::GlobalAssemble(*matrix);


    // Check local
    TrilinosCPPTestExperimentalUtilities::CheckSparseMatrixFromLocalMatrix(*matrix, local_matrix);
}

} // namespace Kratos::Testing
