//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || "__"| | | | "_ \ / _ \/ __|
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


KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalCPPTestUtilitiesGenerateSparseMatrix, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 2 * r_comm.Size();
    auto matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(r_comm, size);
    auto local_matrix = TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalMatrix(size);

    // Check
    TrilinosCPPTestExperimentalUtilities::CheckSparseMatrixFromLocalMatrix(*matrix, local_matrix);
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalCPPTestUtilitiesGenerateSparseMatrixRealCase, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 6;
    std::vector<int> row_indexes = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5};
    std::vector<int> column_indexes = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 2, 3, 4, 5, 2, 3, 4, 5};
    std::vector<double> values = {2069000000.0, 0.0, -2069000000.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2069000000.0, 0.0, 4138000000.0, 0.0, -2069000000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2069000000.0, 0.0, 2069000000.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    auto matrix = TrilinosCPPTestExperimentalUtilities::GenerateSparseMatrix(r_comm, size, row_indexes, column_indexes, values);

    // Check
    TrilinosCPPTestExperimentalUtilities::CheckSparseMatrix(*matrix, row_indexes, column_indexes, values);
}

} // namespace Kratos::Testing
