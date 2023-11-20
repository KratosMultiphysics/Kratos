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

namespace Kratos::Testing
{

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosCPPTestUtilitiesGenerateSparseMatrix, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // Size of the system
    const int size = 6;

    // // Debug info
    // auto matrix = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 0.0, true);
    // std::vector<int> row_indexes;
    // std::vector<int> column_indexes;
    // std::vector<double> values;
    // TrilinosCPPTestUtilities::GenerateSparseMatrixIndexAndValuesVectors(matrix, row_indexes, column_indexes, values, true);

    // Generate matrix
    std::vector<int> row_indexes = {0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5};
    std::vector<int> column_indexes = {0, 1, 0, 1, 2, 1, 2, 3, 2, 3, 4, 3, 4, 5, 4, 5};
    std::vector<double> values = {0.0, -1.0, -1.0, 1.0, -1.0, -1.0, 2.0, -1.0, -1.0, 3.0, -1.0, -1.0, 4.0, -1.0, -1.0, 5.0};
    auto A = TrilinosCPPTestUtilities::GenerateSparseMatrix(r_comm, size, row_indexes, column_indexes, values);

    // Check
    TrilinosCPPTestUtilities::CheckSparseMatrix(A, row_indexes, column_indexes, values);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosCPPTestUtilitiesGenerateSparseMatrixRealCase, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // Size of the system
    const int size = 6;

    // Generate A matrix
    std::vector<int> row_indexes = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5};
    std::vector<int> column_indexes = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 2, 3, 4, 5, 2, 3, 4, 5};
    std::vector<double> values = {2069000000.0, 0.0, -2069000000.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2069000000.0, 0.0, 4138000000.0, 0.0, -2069000000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2069000000.0, 0.0, 2069000000.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    auto A = TrilinosCPPTestUtilities::GenerateSparseMatrix(r_comm, size, row_indexes, column_indexes, values);

    // Check
    TrilinosCPPTestUtilities::CheckSparseMatrix(A, row_indexes, column_indexes, values);

    // Generate T matrix
    row_indexes = {0,1,2,3,4,4,5};
    column_indexes = {0,1,2,3,2,4,5};
    values = {1,1,1,1,1,0,1};
    auto T = TrilinosCPPTestUtilities::GenerateSparseMatrix(r_comm, size, row_indexes, column_indexes, values);

    // Check
    TrilinosCPPTestUtilities::CheckSparseMatrix(T, row_indexes, column_indexes, values);
}

} // namespace Kratos::Testing