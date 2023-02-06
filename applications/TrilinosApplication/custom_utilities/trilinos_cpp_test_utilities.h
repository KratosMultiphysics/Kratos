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

#pragma once

// System includes

// External includes

// Project includes
#include "trilinos_space.h"
#include "includes/checks.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @class TrilinosCPPTestUtilities
 * @brief Utilities to develop C++ tests in Trilinos
 * @ingroup TrilinosApplication
 * @author Vicente Mataix Ferandiz
 */
class TrilinosCPPTestUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TrilinosDofUpdater
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosCPPTestUtilities);

    /// Basic definitions
    using TrilinosSparseSpaceType = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>;
    using TrilinosLocalSpaceType = UblasSpace<double, Matrix, Vector>;

    using TrilinosSparseMatrixType = TrilinosSparseSpaceType::MatrixType;
    using TrilinosVectorType = TrilinosSparseSpaceType::VectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
    * @brief Generates a dummy diagonal sparse matrix for Trilinos
    * @param rDataCommunicator The data communicator considered
    * @param NumGlobalElements The global dimension of the matrix
    * @param Offset The offset considered
    */
    static TrilinosSparseMatrixType GenerateDummySparseMatrix(
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
        TrilinosSparseMatrixType A(Copy, Map, NumNz.data());

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

    /**
    * @brief Generates a dummy vector for Trilinos
    * @param rDataCommunicator The data communicator considered
    * @param NumGlobalElements The global dimension of the matrix
    * @param Offset The offset considered
    */
    static TrilinosVectorType GenerateDummyVector(
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

        // Create a Epetra_Matrix
        TrilinosVectorType b(Map);

        double value;
        for( int i=0 ; i<NumMyElements; ++i ) {
            value = Offset + static_cast<double>(MyGlobalElements[i]);
            b[0][i]= value;
        }

        b.GlobalAssemble();

        return b;
    }

    /**
    * @brief This method checks the values of a sparse matrix with the given indices and values
    * @param rA The matrix to check
    * @param rRowIndexes The row indices
    * @param rColumnIndexes The column indices
    * @param rValues The values
    * @param Tolerance The tolerance considered
    */
    static void CheckSparseMatrix(
        const TrilinosSparseSpaceType::MatrixType& rA,
        const std::vector<int>& rRowIndexes,
        const std::vector<int>& rColumnIndexes,
        const std::vector<double>& rValues,
        const double Tolerance = 1e-8
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

    ///@}

}; /// class TrilinosCPPTestUtilities

} /// namespace Kratos