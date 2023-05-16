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

    using TrilinosLocalMatrixType = TrilinosLocalSpaceType::MatrixType;
    using TrilinosLocalVectorType = TrilinosLocalSpaceType::VectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
    * @brief Generates a dummy diagonal local matrix for Trilinos
    * @param rDataCommunicator The data communicator considered
    * @param NumGlobalElements The global dimension of the matrix
    * @param Offset The offset considered
    * @param AddNoDiagonalValues If adding non diagonal values
    */
    static TrilinosLocalMatrixType GenerateDummyLocalMatrix(
        const int NumGlobalElements = 12,
        const double Offset = 0.0,
        const bool AddNoDiagonalValues = false
        );

    /**
    * @brief Generates a dummy diagonal sparse matrix for Trilinos
    * @param rDataCommunicator The data communicator considered
    * @param NumGlobalElements The global dimension of the matrix
    * @param Offset The offset considered
    * @param AddNoDiagonalValues If adding non diagonal values
    */
    static TrilinosSparseMatrixType GenerateDummySparseMatrix(
        const DataCommunicator& rDataCommunicator,
        const int NumGlobalElements = 12,
        const double Offset = 0.0,
        const bool AddNoDiagonalValues = false
        );

    /**
    * @brief Generates a dummy local vector for Trilinos
    * @param rDataCommunicator The data communicator considered
    * @param NumGlobalElements The global dimension of the matrix
    * @param Offset The offset considered
    */
    static TrilinosLocalVectorType GenerateDummyLocalVector(
        const int NumGlobalElements = 12,
        const double Offset = 0.0
        );

    /**
    * @brief Generates a dummy sparse vector for Trilinos
    * @param rDataCommunicator The data communicator considered
    * @param NumGlobalElements The global dimension of the matrix
    * @param Offset The offset considered
    */
    static TrilinosVectorType GenerateDummySparseVector(
        const DataCommunicator& rDataCommunicator,
        const int NumGlobalElements = 12,
        const double Offset = 0.0
        );

    /**
    * @brief This method checks the values of a sparse vector with the given serial vector
    * @param vector The matrix to check
    * @param vector The reference matrix
    * @param NegligibleValueThreshold The tolerance considered
    */
    static void CheckSparseVectorFromLocalVector(
        const TrilinosVectorType& rA,
        const TrilinosLocalVectorType& rB,
        const double NegligibleValueThreshold = 1e-8
        );

    /**
    * @brief This method checks the values of a sparse vector with the given indices and values
    * @param rb The vector to check
    * @param rIndexes The indices
    * @param rValues The values
    * @param NegligibleValueThreshold The tolerance considered
    */
    static void CheckSparseVector(
        const TrilinosVectorType& rb,
        const std::vector<int>& rIndexes,
        const std::vector<double>& rValues,
        const double NegligibleValueThreshold = 1e-8
        );

    /**
    * @brief This method checks the values of a sparse matrix with the given serial matrix
    * @param rA The matrix to check
    * @param rB The reference matrix
    * @param NegligibleValueThreshold The tolerance considered
    */
    static void CheckSparseMatrixFromLocalMatrix(
        const TrilinosSparseMatrixType& rA,
        const TrilinosLocalMatrixType& rB,
        const double NegligibleValueThreshold = 1e-8
        );

    /**
    * @brief This method checks the values of a sparse matrix with the given indices and values
    * @param rA The matrix to check
    * @param rRowIndexes The row indices
    * @param rColumnIndexes The column indices
    * @param rB The reference matrix
    * @param NegligibleValueThreshold The tolerance considered
    */
    static void CheckSparseMatrixFromLocalMatrix(
        const TrilinosSparseMatrixType& rA,
        const std::vector<int>& rRowIndexes,
        const std::vector<int>& rColumnIndexes,
        const TrilinosLocalMatrixType& rB,
        const double NegligibleValueThreshold = 1e-8
        );

    /**
    * @brief This method checks the values of a sparse matrix with the given indices and values
    * @param rA The matrix to check
    * @param rRowIndexes The row indices
    * @param rColumnIndexes The column indices
    * @param rValues The values
    * @param NegligibleValueThreshold The tolerance considered
    */
    static void CheckSparseMatrix(
        const TrilinosSparseMatrixType& rA,
        const std::vector<int>& rRowIndexes,
        const std::vector<int>& rColumnIndexes,
        const std::vector<double>& rValues,
        const double NegligibleValueThreshold = 1e-8
        );

    /**
    * @brief This method generates a set of row, columns and values from a sparse matrix
    * @param rA The matrix where generate the vectors of rows, indexes and values
    * @param rRowIndexes The row indices
    * @param rColumnIndexes The column indices
    * @param rValues The values
    * @param PrintValues If printing the vectors, for debugging pourposes
    * @param IncludeHardZeros If including the hard zeros (giving  a cerain threshold)
    */
    static void GenerateSparseMatrixIndexAndValuesVectors(
        const TrilinosSparseSpaceType::MatrixType& rA,
        std::vector<int>& rRowIndexes,
        std::vector<int>& rColumnIndexes,
        std::vector<double>& rValues,
        const bool PrintValues = false,
        const double ThresholdIncludeHardZeros = -1
        );

    /**
    * @brief This method generates a sparse matrix from a set of row, columns and values
    * @param rDataCommunicator The data communicator considered
    * @param NumGlobalElements The global dimension of the matrix
    * @param rRowIndexes The row indices
    * @param rColumnIndexes The column indices
    * @param rValues The values
    * @param pMap Map pointer
    * @return The matrix generated
    */
    static TrilinosSparseMatrixType GenerateSparseMatrix(
        const DataCommunicator& rDataCommunicator,
        const int NumGlobalElements,
        const std::vector<int>& rRowIndexes,
        const std::vector<int>& rColumnIndexes,
        const std::vector<double>& rValues,
        const Epetra_Map* pMap =  nullptr
        );

    ///@}

}; /// class TrilinosCPPTestUtilities

} /// namespace Kratos