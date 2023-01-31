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

namespace Kratos::Testing
{

/// The sparse matrix type
using TrilinosSparseSpaceType = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>;
using SparseMatrixType = TrilinosSparseSpaceType::MatrixType;
using VectorType = TrilinosSparseSpaceType::VectorType;

/**
 * @brief Generates a dummy diagonal sparse matrix for Trilinos
 * @param rDataCommunicator The data communicator considered
 * @param NumGlobalElements The global dimension of the matrix
 * @param Offset The offset considered
 */
SparseMatrixType GenerateDummySparseMatrix(
    const DataCommunicator& rDataCommunicator,
    const int NumGlobalElements = 12,
    const double Offset = 0.0
    )
{
    // Generate Epetra communicator
    KRATOS_ERROR_IF_NOT(rDataCommunicator.IsDistributed()) << "Only distributed DataCommunicators can be used!" << std::endl;
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(rDataCommunicator);
    Epetra_MpiComm epetra_comm(raw_mpi_comm);

    // Create a map
    Epetra_Map Map(NumGlobalElements,0,epetra_comm);

    // Local number of rows
    const int NumMyElements = Map.NumMyElements();

    // Get update list
    int* MyGlobalElements = Map.MyGlobalElements( );

    // Create an integer vector NumNz that is used to build the EPetra Matrix.
    // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation
    // on this processor
    std::vector<int> NumNz(NumMyElements, 1);

    // Create a Epetra_Matrix
    SparseMatrixType A(Copy, Map, NumNz.data());

    double value;
    for( int i=0 ; i<NumMyElements; ++i ) {
        // Put in the diagonal entry
        value = Offset + static_cast<double>(MyGlobalElements[i]);
        A.InsertGlobalValues(MyGlobalElements[i], 1, &value, MyGlobalElements+i);
    }

    // Finish up, trasforming the matrix entries into local numbering,
    // to optimize data transfert during matrix-vector products
    A.FillComplete();

    return A;
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosCheckAndCorrectZeroDiagonalValues, KratosTrilinosApplicationMPITestSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");
    auto& r_process_info = r_model_part.GetProcessInfo();
    r_process_info.SetValue(BUILD_SCALE_FACTOR, 1.0);

    // The data communicator
    const DataCommunicator& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 12;
    SparseMatrixType matrix12x12 = GenerateDummySparseMatrix(r_comm, size);

    // Generate Epetra communicator
    KRATOS_ERROR_IF_NOT(r_comm.IsDistributed()) << "Only distributed DataCommunicators can be used!" << std::endl;
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    Epetra_MpiComm epetra_comm(raw_mpi_comm);

    // Dummy vector
    Epetra_Map map(size, 0, epetra_comm);
    VectorType vector12(map);

    // Test the norm of the matrix
    double norm = 0.0;
    norm = TrilinosSparseSpaceType::CheckAndCorrectZeroDiagonalValues(r_process_info, matrix12x12, vector12, SCALING_DIAGONAL::NO_SCALING, r_comm);
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
    const DataCommunicator& r_comm = Testing::GetDefaultDataCommunicator();

    // The dummy matrix
    const int size = 12;
    SparseMatrixType matrix12x12 = GenerateDummySparseMatrix(r_comm, size, 1.0);

    // Test the norm of the matrix
    double norm = 0.0;
    norm = TrilinosSparseSpaceType::GetScaleNorm(r_process_info, matrix12x12, SCALING_DIAGONAL::NO_SCALING, r_comm);
    KRATOS_CHECK_DOUBLE_EQUAL(norm, 1.0);
    norm = TrilinosSparseSpaceType::GetScaleNorm(r_process_info, matrix12x12, SCALING_DIAGONAL::CONSIDER_PRESCRIBED_DIAGONAL, r_comm);
    KRATOS_CHECK_DOUBLE_EQUAL(norm, 3.0);
    norm = TrilinosSparseSpaceType::GetScaleNorm(r_process_info, matrix12x12, SCALING_DIAGONAL::CONSIDER_NORM_DIAGONAL, r_comm);
    KRATOS_CHECK_NEAR(norm, 2.124591464, 1.0e-6);
    norm = TrilinosSparseSpaceType::GetScaleNorm(r_process_info, matrix12x12, SCALING_DIAGONAL::CONSIDER_MAX_DIAGONAL, r_comm);
    KRATOS_CHECK_DOUBLE_EQUAL(norm, 12.0);
    norm = TrilinosSparseSpaceType::GetAveragevalueDiagonal(matrix12x12, r_comm);
    KRATOS_CHECK_DOUBLE_EQUAL(norm, 6.5);
    norm = TrilinosSparseSpaceType::GetMinDiagonal(matrix12x12, r_comm);
    KRATOS_CHECK_DOUBLE_EQUAL(norm, 1.0);
}

} // namespace Kratos::Testing
