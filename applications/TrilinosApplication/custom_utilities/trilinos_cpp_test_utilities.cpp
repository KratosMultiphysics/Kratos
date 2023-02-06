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
#include "includes/checks.h"
#include "custom_utilities/trilinos_cpp_test_utilities.h"

namespace Kratos
{

TrilinosCPPTestUtilities::TrilinosLocalMatrixType TrilinosCPPTestUtilities::GenerateDummyLocalMatrix(
    const int NumGlobalElements,
    const double Offset,
    const bool AddNoDiagonalValues
    )
{
    TrilinosLocalMatrixType A(NumGlobalElements, NumGlobalElements);
    A.clear();

    for (int i = 0; i < NumGlobalElements; ++i) {
        A(i, i) = Offset + static_cast<double>(i);
    }
    if (AddNoDiagonalValues) {
        for (int i = 0; i < NumGlobalElements; ++i) {
            if (i - 1 >= 0) A(i - 1, i) = -1;
            if (i + 1 < NumGlobalElements) A(i + 1, i) = -1;
        }
    }

    return A;
}

/***********************************************************************************/
/***********************************************************************************/

TrilinosCPPTestUtilities::TrilinosSparseMatrixType TrilinosCPPTestUtilities::GenerateDummySparseMatrix(
    const DataCommunicator& rDataCommunicator,
    const int NumGlobalElements,
    const double Offset,
    const bool AddNoDiagonalValues
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
    if (AddNoDiagonalValues) {
        // So we need 2 off-diagonal terms (except for the first and last equation)
        for ( int i=0; i<NumMyElements; i++) {
            if (MyGlobalElements[i]==0 || MyGlobalElements[i] == NumGlobalElements-1) {
                NumNz[i] = 2;
            } else {
                NumNz[i] = 3;
            }
        }
    }

    // Create an Epetra_Matrix
    TrilinosSparseMatrixType A(Copy, Map, NumNz.data());

    double* non_diagonal_values = new double[2];
    non_diagonal_values[0] = -1.0; non_diagonal_values[1] = -1.0;
    int* non_diagonal_indices = new int[2];
    int num_non_diagonal_entries;
    double value;
    for( int i=0 ; i<NumMyElements; ++i ) {
        // Non diagonal values
        if (AddNoDiagonalValues) {
            if (MyGlobalElements[i]==0) {
                non_diagonal_indices[0] = 1;
                num_non_diagonal_entries = 1;
            } else if (MyGlobalElements[i] == NumGlobalElements-1) {
                non_diagonal_indices[0] = NumGlobalElements-2;
                num_non_diagonal_entries = 1;
            } else {
                non_diagonal_indices[0] = MyGlobalElements[i]-1;
                non_diagonal_indices[1] = MyGlobalElements[i]+1;
                num_non_diagonal_entries = 2;
            }
            A.InsertGlobalValues(MyGlobalElements[i], num_non_diagonal_entries, non_diagonal_values, non_diagonal_indices);
        }
        // Put in the diagonal entry
        value = Offset + static_cast<double>(MyGlobalElements[i]);
        A.InsertGlobalValues(MyGlobalElements[i], 1, &value, MyGlobalElements+i);
    }

    // Finish up, trasforming the matrix entries into local numbering,
    // to optimize data transfert during matrix-vector products
    A.FillComplete();

    return A;
}

/***********************************************************************************/
/***********************************************************************************/

TrilinosCPPTestUtilities::TrilinosLocalVectorType TrilinosCPPTestUtilities::GenerateDummyLocalVector(
    const int NumGlobalElements,
    const double Offset
    )
{
    TrilinosLocalVectorType b(NumGlobalElements);
    for (int i = 0; i < NumGlobalElements; ++i) {
        b[i] = Offset + static_cast<double>(i);
    }

    return b;
}

/***********************************************************************************/
/***********************************************************************************/

TrilinosCPPTestUtilities::TrilinosVectorType TrilinosCPPTestUtilities::GenerateDummySparseVector(
    const DataCommunicator& rDataCommunicator,
    const int NumGlobalElements,
    const double Offset
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

    // Create an Epetra_Vector
    TrilinosVectorType b(Map);

    double value;
    for( int i=0 ; i<NumMyElements; ++i ) {
        value = Offset + static_cast<double>(MyGlobalElements[i]);
        b[0][i]= value;
    }

    b.GlobalAssemble();

    return b;
}

/***********************************************************************************/
/***********************************************************************************/

void TrilinosCPPTestUtilities::CheckSparseVectorFromLocalVector(
    const TrilinosVectorType& rA,
    const TrilinosLocalVectorType& rB,
    const double Tolerance
    )
{
    const std::size_t total_size = rB.size();
    std::vector<int> indexes;
    indexes.reserve(total_size);
    std::vector<double> values;
    values.reserve(total_size);
    for (std::size_t i = 0; i < rB.size(); ++i) {
        indexes.push_back(i);
        values.push_back(rB[i]);
    }
    CheckSparseVector(rA, indexes, values);
}

/***********************************************************************************/
/***********************************************************************************/

void TrilinosCPPTestUtilities::CheckSparseVector(
    const TrilinosVectorType& rb,
    const std::vector<int>& rIndexes,
    const std::vector<double>& rValues,
    const double Tolerance
    )
{
    int index;
    double value;
    const auto& r_map = rb.Map();
    for (std::size_t counter = 0; counter < rIndexes.size(); ++counter) {
        index = rIndexes[counter];
        value = rValues[counter];
        if (r_map.MyGID(index)) {
            const double ref_value = rb[0][r_map.LID(index)];
            KRATOS_CHECK_RELATIVE_NEAR(value, ref_value, Tolerance);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void TrilinosCPPTestUtilities::CheckSparseMatrixFromLocalMatrix(
    const TrilinosSparseMatrixType& rA,
    const TrilinosLocalMatrixType& rB,
    const double Tolerance
    )
{
    const std::size_t total_size = rB.size1() * rB.size2();
    std::vector<int> row_indexes;
    row_indexes.reserve(total_size);
    std::vector<int> column_indexes;
    column_indexes.reserve(total_size);
    std::vector<double> values;
    values.reserve(total_size);
    for (std::size_t i = 0; i < rB.size1(); ++i) {
        for (std::size_t j = 0; j < rB.size2(); ++j) {
            row_indexes.push_back(i);
            column_indexes.push_back(j);
            values.push_back(rB(i, j));
        }
    }
    CheckSparseMatrix(rA, row_indexes, column_indexes, values);
}

/***********************************************************************************/
/***********************************************************************************/

void TrilinosCPPTestUtilities::CheckSparseMatrix(
    const TrilinosSparseMatrixType& rA,
    const std::vector<int>& rRowIndexes,
    const std::vector<int>& rColumnIndexes,
    const std::vector<double>& rValues,
    const double Tolerance
    )
{
    int row, column;
    double value;
    for (std::size_t counter = 0; counter < rRowIndexes.size(); ++counter) {
        row = rRowIndexes[counter];
        column = rColumnIndexes[counter];
        value = rValues[counter];
        for (int i = 0; i < rA.NumMyRows(); i++) {
            int numEntries; // Number of non-zero entries
            double* vals;   // Row non-zero values
            int* cols;      // Column indices of row non-zero values
            rA.ExtractMyRowView(i, numEntries, vals, cols);
            const int row_gid = rA.RowMap().GID(i);
            if (row == row_gid) {
                int j;
                for (j = 0; j < numEntries; j++) {
                    const int col_gid = rA.ColMap().GID(cols[j]);
                    if (col_gid == column) {
                        KRATOS_CHECK_RELATIVE_NEAR(value, vals[j], Tolerance)
                        break;
                    }
                }
                break;
            }
        }
    }
}

} /// namespace Kratos