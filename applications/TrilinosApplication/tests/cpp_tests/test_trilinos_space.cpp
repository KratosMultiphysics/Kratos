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
#include "custom_utilities/trilinos_cpp_test_utilities.h"

namespace Kratos::Testing
{

/// The sparse matrix type
using TrilinosSparseSpaceType = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>;
using SparseMatrixType = TrilinosSparseSpaceType::MatrixType;
using VectorType = TrilinosSparseSpaceType::VectorType;

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosSizeVector, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto vector = TrilinosCPPTestUtilities::GenerateDummyVector(r_comm, size);

    // Check
    KRATOS_CHECK_EQUAL(size, TrilinosSparseSpaceType::Size(vector));
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosSizeMatrix, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto vector = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size);

    // Check
    KRATOS_CHECK_EQUAL(size, TrilinosSparseSpaceType::Size1(vector));
    KRATOS_CHECK_EQUAL(size, TrilinosSparseSpaceType::Size2(vector));
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosDotProduct, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto vector1 = TrilinosCPPTestUtilities::GenerateDummyVector(r_comm, size);
    auto vector2 = TrilinosCPPTestUtilities::GenerateDummyVector(r_comm, size);
    double ref = 0.0;
    for (int i = 0; i < size; ++i) {
        ref += std::pow(static_cast<double>(i), 2);
    }

    // Check
    KRATOS_CHECK_DOUBLE_EQUAL(ref, TrilinosSparseSpaceType::Dot(vector1, vector2));
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosMaxMin, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto vector = TrilinosCPPTestUtilities::GenerateDummyVector(r_comm, size);

    // Check
    KRATOS_CHECK_DOUBLE_EQUAL(0.0, TrilinosSparseSpaceType::Min(vector));
    KRATOS_CHECK_DOUBLE_EQUAL(static_cast<double>(size - 1), TrilinosSparseSpaceType::Max(vector));
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosTwoNormVector, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy vector
    const int size = 2 * r_comm.Size();
    auto vector = TrilinosCPPTestUtilities::GenerateDummyVector(r_comm, size);
    double ref = 0.0;
    for (int i = 0; i < size; ++i) {
        ref += std::pow(static_cast<double>(i), 2);
    }
    ref = std::sqrt(ref);

    // Check
    KRATOS_CHECK_DOUBLE_EQUAL(ref, TrilinosSparseSpaceType::TwoNorm(vector));
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosTwoNormMatrix, KratosTrilinosApplicationMPITestSuite)
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
    KRATOS_CHECK_DOUBLE_EQUAL(ref, TrilinosSparseSpaceType::TwoNorm(matrix));
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
    VectorType vector12(map);

    // Test the norm of the matrix
    double norm = 0.0;
    norm = TrilinosSparseSpaceType::CheckAndCorrectZeroDiagonalValues(r_process_info, matrix12x12, vector12, r_comm, SCALING_DIAGONAL::NO_SCALING);
    KRATOS_CHECK_DOUBLE_EQUAL(norm, 1.0);
    if (r_comm.Rank() == 0) {
        const auto& raw_values = matrix12x12.ExpertExtractValues();
        KRATOS_CHECK_DOUBLE_EQUAL(raw_values[0], 1.0);
    }
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
    norm = TrilinosSparseSpaceType::GetScaleNorm(r_process_info, matrix12x12, r_comm, SCALING_DIAGONAL::NO_SCALING);
    KRATOS_CHECK_DOUBLE_EQUAL(norm, 1.0);
    norm = TrilinosSparseSpaceType::GetScaleNorm(r_process_info, matrix12x12, r_comm, SCALING_DIAGONAL::CONSIDER_PRESCRIBED_DIAGONAL);
    KRATOS_CHECK_DOUBLE_EQUAL(norm, 3.0);
    norm = TrilinosSparseSpaceType::GetScaleNorm(r_process_info, matrix12x12, r_comm, SCALING_DIAGONAL::CONSIDER_NORM_DIAGONAL);
    KRATOS_CHECK_NEAR(norm, 2.124591464, 1.0e-6);
    norm = TrilinosSparseSpaceType::GetScaleNorm(r_process_info, matrix12x12, r_comm, SCALING_DIAGONAL::CONSIDER_MAX_DIAGONAL);
    KRATOS_CHECK_DOUBLE_EQUAL(norm, 12.0);
    norm = TrilinosSparseSpaceType::GetAveragevalueDiagonal(matrix12x12, r_comm);
    KRATOS_CHECK_DOUBLE_EQUAL(norm, 6.5);
    norm = TrilinosSparseSpaceType::GetMinDiagonal(matrix12x12, r_comm);
    KRATOS_CHECK_DOUBLE_EQUAL(norm, 1.0);
}

} // namespace Kratos::Testing