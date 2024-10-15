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

/* Trilinos includes */
#include <Tpetra_Core.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_FECrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <TpetraExt_MatrixMatrix.hpp>
// #include <TpetraExt_TripleMatrixMultiply_def.hpp> # TODO: Failing compilation for some reason in Ubuntu 20.04 version available
//#include <MatrixMarket_Tpetra.hpp> # TODO: Failing compilation for some reason in Ubuntu 20.04 version available

// Project includes
#include "includes/ublas_interface.h"
#include "spaces/ublas_space.h"
#include "includes/data_communicator.h"
#include "mpi/includes/mpi_data_communicator.h"
//#include "custom_utilities/trilinos_dof_updater.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class TrilinosSpaceExperimental
 * @ingroup TrilinosApplication
 * @brief The space adapted for Trilinos vectors and matrices (TPetra)
 * @details This class is experimental and aims to bring support for new TPetra matrices and vectors to TrilinosApplication
 * @author Vicente Mataix Ferrandiz
 * @tparam TMatrixType The matrix type considered
 * @tparam TVectorType the vector type considered
 */
template<class TMatrixType, class TVectorType>
class TrilinosSpaceExperimental
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TrilinosSpaceExperimental
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosSpaceExperimental);

    /// Definition of the data type
    using DataType = typename TMatrixType::scalar_type;

    /// Definition of the matrix type
    using MatrixType = TMatrixType;

    /// Definition of the vector type
    using VectorType = TVectorType;

    /// Definition of the index type
    using IndexType = std::size_t;

    /// Definition of the size type
    using SizeType = std::size_t;

    /// Class definition
    using ClassType = TrilinosSpaceExperimental<TMatrixType, TVectorType>;

    /// Tpetra definitions
    // Your scalar type; the type of sparse matrix entries. e.g., double.
    using ST = typename MatrixType::scalar_type;
    // Your local ordinal type; the signed integer type used to store local sparse matrix indices.  e.g., int.
    using LO = typename MatrixType::local_ordinal_type;
    // Your global ordinal type; the signed integer type used to index the matrix globally, over all processes. e.g., int, long, ptrdif_t, int64_t, ...
    using GO = typename MatrixType::global_ordinal_type;
    // The Node type.  e.g., Kokkos::DefaultNode::DefaultNodeType, defined in KokkosCompat_DefaultNode.hpp.
    using NT = typename MatrixType::node_type;

    /// Define the import/export types
    using ImportType = Tpetra::Import<LO, GO, NT>;
    using ExportType = Tpetra::Export<LO, GO, NT>;

    /// Define the map type
    using MapType = Tpetra::Map<LO, GO, NT>;
    using MapPointerType = Teuchos::RCP<const MapType>;

    /// Define the graph type
    using GraphType = Tpetra::FECrsGraph<LO, GO, NT>;
    using GraphPointerType = Teuchos::RCP<const GraphType>;

    // Define TPetra communicator
    using CommunicatorType = Teuchos::MpiComm<int>;
    using CommunicatorPointerType = Teuchos::RCP<const CommunicatorType>;

    /// Definition of the pointer types
    using MatrixPointerType = Teuchos::RCP<MatrixType>;
    using VectorPointerType = Teuchos::RCP<VectorType>;

    /// Column indices of row non-zero values type
    using ColumnViewType = typename MatrixType::nonconst_global_inds_host_view_type;

    /// Row non-zero values type
    using ValueViewType = typename MatrixType::nonconst_values_host_view_type;

    /// Some other definitions
    // using DofUpdaterType = TrilinosDofUpdater<ClassType>;
    // using DofUpdaterPointerType = typename DofUpdater<ClassType>::UniquePointer;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    TrilinosSpaceExperimental()
    {
    }

    /// Destructor.
    virtual ~TrilinosSpaceExperimental()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method creates an empty pointer to a matrix
     * @return The pointer to the matrix
     */
    static MatrixPointerType CreateEmptyMatrixPointer()
    {
        return MatrixPointerType(nullptr);
    }

    /**
     * @brief This method creates an empty pointer to a vector
     * @return The pointer to the vector
     */
    static VectorPointerType CreateEmptyVectorPointer()
    {
        return VectorPointerType(nullptr);
    }

    /**
     * @brief This method creates an empty pointer to a matrix using TPetra communicator
     * @param rComm The Tpetra communicator
     * @return The pointer to the matrix
     */
    static MatrixPointerType CreateEmptyMatrixPointer(CommunicatorPointerType pComm)
    {
        const int global_elems = 0;
        MapPointerType map = Teuchos::rcp(new MapType(global_elems, 0, pComm));
        GraphPointerType graph = Teuchos::rcp(new GraphType(map, map, 0));
        return Teuchos::rcp(new MatrixType(graph));
    }

    /**
     * @brief This method creates an empty pointer to a vector using Tpetra communicator
     * @param pComm The Tpetra communicator
     * @return The pointer to the vector
     */
    static VectorPointerType CreateEmptyVectorPointer(CommunicatorPointerType pComm)
    {
        const int global_elems = 0;
        MapPointerType map = Teuchos::rcp(new MapType(global_elems, 0, pComm));
        return Teuchos::rcp(new VectorType(map));
    }

    /**
     * @brief Returns size of vector rV
     * @param rV The vector considered
     * @return The size of the vector
     */
    static IndexType Size(const VectorType& rV)
    {
        const int size = rV.getGlobalLength();
        return size;
    }

    /**
     * @brief Returns number of rows of rM
     * @param rM The matrix considered
     * @return The number of rows of rM
     */
    static IndexType Size1(MatrixType const& rM)
    {
        const int size1 = rM.getGlobalNumRows();
        return size1;
    }

    /**
     * @brief Returns number of columns of rM
     * @param rM The matrix considered
     * @return The number of columns of rM
     */
    static IndexType Size2(MatrixType const& rM)
    {
        const int size1 = rM.getGlobalNumCols();
        return size1;
    }

    /**
     * @brief Returns the column of the matrix in the given position
     * @details rXi = rMij
     * @param j The position of the column
     * @param rM The matrix considered
     * @param rX The column considered
     * @todo Implement this method
     */
    static void GetColumn(
        const unsigned int j,
        const MatrixType& rM,
        VectorType& rX
        )
    {
        KRATOS_ERROR << "GetColumn method is not currently implemented" << std::endl;
    }

    /**
     * @brief Returns a copy of the matrix rX
     * @details rY = rX
     * @param rX The matrix considered
     * @param rY The copy of the matrix rX
     */
    static void Copy(
        const MatrixType& rX,
        MatrixType& rY
        )
    {
        rY = rX;
    }

    /**
     * @brief Returns a copy of the vector rX
     * @details rY = rX
     * @param rX The vector considered
     * @param rY The copy of the vector rX
     */
    static void Copy(
        const VectorType& rX,
        VectorType& rY
        )
    {
        rY = rX;
    }

    /**
     * @brief Returns the product of two vectors
     * @details rX * rY
     * @param rX The first vector considered
     * @param rY The second vector considered
     */
    static double Dot(
        const VectorType& rX,
        const VectorType& rY
        )
    {
        return rX.dot(rY);
    }

    /**
     * @brief Returns the maximum value of the vector rX
     * @param rX The vector considered
     * @return The maximum value of the vector rX
     */
    static double Max(const VectorType& rX)
    {
        // Access the local data
        auto localVec = rX.getLocalViewHost(Tpetra::Access::ReadOnly);

        // Find the local maximum
        ST localMax = localVec(0,0);
        auto localLength = rX.getLocalLength();
        for (std::size_t i = 1; i < localLength; ++i) {
            localMax = std::max(localMax, localVec(i,0));
        }

        // Perform a global maximum reduction
        ST globalMax;
        Teuchos::reduceAll(*rX.getMap()->getComm(), Teuchos::REDUCE_MAX, localMax, Teuchos::outArg(globalMax));

        return globalMax;
    }

    /**
     * @brief Returns the minimum value of the vector rX
     * @param rX The vector considered
     * @return The minimum value of the vector rX
     */
    static double Min(const VectorType& rX)
    {
        // Access the local data
        auto localVec = rX.getLocalViewHost(Tpetra::Access::ReadOnly);

        // Find the local minimum
        ST localMin = localVec(0,0);
        auto localLength = rX.getLocalLength();
        for (std::size_t i = 1; i < localLength; ++i) {
            localMin = std::min(localMin, localVec(i,0));
        }

        // Perform a global minimum reduction
        ST globalMin;
        Teuchos::reduceAll(*rX.getMap()->getComm(), Teuchos::REDUCE_MIN, localMin, Teuchos::outArg(globalMin));

        return globalMin;
    }

    /**
     * @brief Returns the norm of the vector rX
     * @details ||rX||2
     * @param rX The vector considered
     * @return The norm of the vector rX
     */
    static double TwoNorm(const VectorType& rX)
    {
        // Use the built-in norm2() method to compute the Euclidean norm (2-norm)
        return rX.norm2();
    }

    /**
     * @brief Returns the Frobenius norm of the matrix rA
     * @details ||rA||2
     * @param rA The matrix considered
     * @return The Frobenius norm of the matrix rA
     */
    static double TwoNorm(const MatrixType& rA)
    {
        // Use the built-in getFrobeniusNorm() method to compute the Frobenius norm
        return rA.getFrobeniusNorm();
    }

    /**
     * @brief Returns the multiplication of a matrix by a vector
     * @details y = A*x
     * @param rA The matrix considered
     * @param rX The vector considered
     * @param rY The result of the multiplication
     */
    static void Mult(
        const MatrixType& rA,
        const VectorType& rX,
        VectorType& rY
        )
    {
        // Multiply: y = A * x
        rA.apply(rX, rY);
    }

    /**
     * @brief Returns the multiplication matrix-matrix
     * @details C = A*B
     * @param rA The first matrix considered
     * @param rB The second matrix considered
     * @param rC The result of the multiplication
     * @param CallFillCompleteOnResult	Optional argument, defaults to true. Power users may specify this argument to be false if they DON'T want this function to call C.FillComplete. (It is often useful to allow this function to call C.FillComplete, in cases where one or both of the input matrices are rectangular and it is not trivial to know which maps to use for the domain- and range-maps.)
     * @param KeepAllHardZeros	Optional argument, defaults to false. If true, Multiply, keeps all entries in C corresponding to hard zeros. If false, the following happens by case: A*B^T, A^T*B^T - Does not store entries caused by hard zeros in C. A^T*B (unoptimized) - Hard zeros are always stored (this option has no effect) A*B, A^T*B (optimized) - Hard zeros in corresponding to hard zeros in A are not stored, There are certain cases involving reuse of C, where this can be useful.
     */
    static void Mult(
        const MatrixType& rA,
        const MatrixType& rB,
        MatrixType& rC,
        const bool CallFillCompleteOnResult = true,
        const bool KeepAllHardZeros = false
        )
    {
        KRATOS_TRY

        // Check if the KeepAllHardZeros is implemented
        //KRATOS_ERROR_IF(KeepAllHardZeros) << "KeepAllHardZeros is not implemented in Trilinos TPetra" << std::endl;

        // Ensure rC is ready to receive the result
        rC.resumeFill();

        // Perform the matrix multiplication
        Tpetra::MatrixMatrix::Multiply(rA, false, rB, false, rC, CallFillCompleteOnResult);

        KRATOS_CATCH("")
    }

    /**
     * @brief Returns the transpose multiplication of a matrix by a vector
     * @details y = AT*x
     * @param rA The matrix considered
     * @param rX The vector considered
     * @param rY The result of the multiplication
     */
    static void TransposeMult(
        const MatrixType& rA,
        const VectorType& rX,
        VectorType& rY
        )
    {
        // Use the apply method with the transpose flag set to true
        rA.apply(rX, rY, Teuchos::TRANS);
    }

    /**
     * @brief Returns the transpose multiplication matrix-matrix
     * @details C = A*B
     * @param rA The first matrix considered
     * @param rB The second matrix considered
     * @param rC The result of the multiplication
     * @param TransposeFlag Flags to transpose the matrices
     * @param CallFillCompleteOnResult	Optional argument, defaults to true. Power users may specify this argument to be false if they DON'T want this function to call C.FillComplete. (It is often useful to allow this function to call C.FillComplete, in cases where one or both of the input matrices are rectangular and it is not trivial to know which maps to use for the domain- and range-maps.)
     * @param KeepAllHardZeros	Optional argument, defaults to false. If true, Multiply, keeps all entries in C corresponding to hard zeros. If false, the following happens by case: A*B^T, A^T*B^T - Does not store entries caused by hard zeros in C. A^T*B (unoptimized) - Hard zeros are always stored (this option has no effect) A*B, A^T*B (optimized) - Hard zeros in corresponding to hard zeros in A are not stored, There are certain cases involving reuse of C, where this can be useful.
     */
    static void TransposeMult(
        const MatrixType& rA,
        const MatrixType& rB,
        MatrixType& rC,
        const std::pair<bool, bool> TransposeFlag = {false, false},
        const bool CallFillCompleteOnResult = true,
        const bool KeepAllHardZeros = false
        )
    {
        KRATOS_TRY

        // Check if the KeepAllHardZeros is implemented
        //KRATOS_ERROR_IF(KeepAllHardZeros) << "KeepAllHardZeros is not implemented in Trilinos TPetra" << std::endl;

        // Ensure rC is ready to receive the result
        rC.resumeFill();

        // Perform the matrix multiplication
        Tpetra::MatrixMatrix::Multiply(rA, TransposeFlag.first, rB, TransposeFlag.second, rC, CallFillCompleteOnResult);

        KRATOS_CATCH("")
    }

    /**
     * @brief Calculates the product operation B'DB
     * @param rA The resulting matrix
     * @param rD The "center" matrix
     * @param rB The matrices to be transposed
     * @param CallFillCompleteOnResult	Optional argument, defaults to true. Power users may specify this argument to be false if they DON'T want this function to call C.FillComplete. (It is often useful to allow this function to call C.FillComplete, in cases where one or both of the input matrices are rectangular and it is not trivial to know which maps to use for the domain- and range-maps.)
     * @param KeepAllHardZeros	Optional argument, defaults to false. If true, Multiply, keeps all entries in C corresponding to hard zeros. If false, the following happens by case: A*B^T, A^T*B^T - Does not store entries caused by hard zeros in C. A^T*B (unoptimized) - Hard zeros are always stored (this option has no effect) A*B, A^T*B (optimized) - Hard zeros in corresponding to hard zeros in A are not stored, There are certain cases involving reuse of C, where this can be useful.
     * @param EnforceInitialGraph If the initial graph is enforced, or a new one is generated
     * @todo TpetraExt_TripleMatrixMultiply_def is failing compilation in the version of Trilinos available in Ubuntu 20.04. We cannot use TripleMatrixMultiply::MultiplyRAP
     */
    static void BtDBProductOperation(
        MatrixType& rA,
        const MatrixType& rD,
        const MatrixType& rB,
        const bool CallFillCompleteOnResult = true,
        const bool KeepAllHardZeros = false,
        const bool EnforceInitialGraph = false
        )
    {
        // Check if the KeepAllHardZeros is implemented
        //KRATOS_ERROR_IF(KeepAllHardZeros) << "KeepAllHardZeros is not implemented in Trilinos TPetra" << std::endl;

        // Check if the EnforceInitialGraph is implemented
        //KRATOS_ERROR_IF(EnforceInitialGraph) << "EnforceInitialGraph is not implemented in Trilinos TPetra" << std::endl;

        // Ensure rA is empty and ready to be filled
        rA.resumeFill();

        // Define first auxiliary matrix
        const auto& r_row_map_a = rA.getRowMap();
        GraphPointerType graph_a = Teuchos::rcp(new GraphType(r_row_map_a, r_row_map_a, 0));
        MatrixPointerType aux_1 = Teuchos::rcp(new MatrixType(graph_a));

        // First multiplication
        TransposeMult(rB, rD, *aux_1, {true, false}, CallFillCompleteOnResult, KeepAllHardZeros);

        // Already existing matrix
        if (rA.getGlobalNumEntries() > 0) {
            // Create a Tpetra_Matrix
            const auto& r_row_map_b = rB.getRowMap();
            GraphPointerType graph_b = Teuchos::rcp(new GraphType(r_row_map_b, r_row_map_b, 0));
            MatrixPointerType aux_2 = Teuchos::rcp(new MatrixType(graph_b));

            // Second multiplication
            Mult(*aux_1, rB, *aux_2, CallFillCompleteOnResult, KeepAllHardZeros);

            // Doing a swap
            MatrixPointerType p_A = Teuchos::rcp(&rA, false);
            std::swap(p_A, aux_2);
        } else { // Empty matrix
            // Second multiplication
            Mult(*aux_1, rB, rA, CallFillCompleteOnResult, KeepAllHardZeros);
        }

    }

    /**
     * @brief Calculates the product operation BDB'
     * @param rA The resulting matrix
     * @param rD The "center" matrix
     * @param rB The matrices to be transposed
     * @param CallFillCompleteOnResult	Optional argument, defaults to true. Power users may specify this argument to be false if they DON'T want this function to call C.FillComplete. (It is often useful to allow this function to call C.FillComplete, in cases where one or both of the input matrices are rectangular and it is not trivial to know which maps to use for the domain- and range-maps.)
     * @param KeepAllHardZeros	Optional argument, defaults to false. If true, Multiply, keeps all entries in C corresponding to hard zeros. If false, the following happens by case: A*B^T, A^T*B^T - Does not store entries caused by hard zeros in C. A^T*B (unoptimized) - Hard zeros are always stored (this option has no effect) A*B, A^T*B (optimized) - Hard zeros in corresponding to hard zeros in A are not stored, There are certain cases involving reuse of C, where this can be useful.
     * @param EnforceInitialGraph If the initial graph is enforced, or a new one is generated
     */
    static void BDBtProductOperation(
        MatrixType& rA,
        const MatrixType& rD,
        const MatrixType& rB,
        const bool CallFillCompleteOnResult = true,
        const bool KeepAllHardZeros = false,
        const bool EnforceInitialGraph = false
        )
    {
        // Check if the KeepAllHardZeros is implemented
        //KRATOS_ERROR_IF(KeepAllHardZeros) << "KeepAllHardZeros is not implemented in Trilinos TPetra" << std::endl;

        // Check if the EnforceInitialGraph is implemented
        //KRATOS_ERROR_IF(EnforceInitialGraph) << "EnforceInitialGraph is not implemented in Trilinos TPetra" << std::endl;

        // Ensure rA is empty and ready to be filled
        rA.resumeFill();

        // Define first auxiliary matrix
        const auto& r_row_map_a = rA.getRowMap();
        GraphPointerType graph_a = Teuchos::rcp(new GraphType(r_row_map_a, r_row_map_a, 0));
        MatrixPointerType aux_1 = Teuchos::rcp(new MatrixType(graph_a));

        // First multiplication
        Mult(rB, rD, *aux_1, CallFillCompleteOnResult, KeepAllHardZeros);

        // Already existing matrix
        if (rA.getGlobalNumEntries() > 0) {
            // Create a Tpetra_Matrix
            const auto& r_row_map_b = rB.getRowMap();
            GraphPointerType graph_b = Teuchos::rcp(new GraphType(r_row_map_b, r_row_map_b, 0));
            MatrixPointerType aux_2 = Teuchos::rcp(new MatrixType(graph_b));

            // Second multiplication
            TransposeMult(*aux_1, rB, *aux_2, {false, true}, CallFillCompleteOnResult, KeepAllHardZeros);

            // Doing a swap
            MatrixPointerType p_A = Teuchos::rcp(&rA, false);
            std::swap(p_A, aux_2);
        } else { // Empty matrix
            // Second multiplication
            TransposeMult(*aux_1, rB, rA, {false, true}, CallFillCompleteOnResult, KeepAllHardZeros);
        }
}

    /**
     * @brief Returns the multiplication of a vector by a scalar
     * @details y = A*x
     * Checks if a multiplication is needed and tries to do otherwise
     * @param rX The vector considered
     * @param A The scalar considered
     */
    static void InplaceMult(
        VectorType& rX,
        const double A
        )
    {
        if (A != 1.00) {
            // Scale the vector: x = A * x
            rX.scale(A);
        }
    }

    /**
     * @brief Returns the multiplication of a vector by a scalar
     * @details x = A*y
     * Checks if a multiplication is needed and tries to do otherwise
     * @note ATTENTION it is assumed no aliasing between rX and rY
     * @param rX The resulting vector considered
     * @param A The scalar considered
     * @param rY The multiplied vector considered
     */
    static void Assign(
        VectorType& rX,
        const double A,
        const VectorType& rY
        )
    {
        if (A != 1.00) {
            // Perform the operation x = A * y
            rX.update(A, rY, 0.0);
        } else {
            rX = rY;
        }
    }

    /**
     * @brief Returns the unaliased addition of a vector by a scalar times a vector
     * @details X += A*y;
     * Checks if a multiplication is needed and tries to do otherwise
     * @note ATTENTION it is assumed no aliasing between rX and rY
     * @param rX The resulting vector considered
     * @param A The scalar considered
     * @param rY The multiplied vector considered
     */
    static void UnaliasedAdd(
        VectorType& rX,
        const double A,
        const VectorType& rY
        )
    {
        rX.update(A, rY, 1.0);
    }

    /**
     * @brief Returns the unaliased addition of two vectors by a scalar
     * @details rZ = (A * rX) + (B * rY)
     * @param A The scalar considered
     * @param rX The first vector considered
     * @param B The scalar considered
     * @param rY The second vector considered
     * @param rZ The resulting vector considered
     */
    static void ScaleAndAdd(
        const double A,
        const VectorType& rX,
        const double B,
        const VectorType& rY,
        VectorType& rZ
        )
    {
        // Compute rZ = A * rX + B * rY
        rZ.update(A, rX, B, rY, 0.0);
    }

    /**
     * @brief Returns the unaliased addition of two vectors by a scalar
     * @details rY = (A * rX) + (B * rY)
     * @param A The scalar considered
     * @param rX The first vector considered
     * @param B The scalar considered
     * @param rY The resulting vector considered
     */
    static void ScaleAndAdd(
        const double A,
        const VectorType& rX,
        const double B,
        VectorType& rY
        )
    {
        // Compute rY = A * rX + B * rY
        rY.update(A, rX, B);
    }

    /**
     * @brief Sets a value in a vector
     * @param rX The vector considered
     * @param i The index of the value considered
     * @param value The value considered
     */
    static void SetValue(
        VectorType& rX,
        IndexType i,
        const double value
        )
    {
        // Get the local index corresponding to the global index `i`
        auto map = rX.getMap();
        IndexType localIndex = map->getLocalElement(i);

        // Check if the index is on this process
        if (localIndex != Tpetra::Details::OrdinalTraits<IndexType>::invalid()) {
            // Set the value at the specified local index
            rX.replaceLocalValue(localIndex, value);
        }
        // If the index `i` is not local, it is ignored (handled by Tpetra's parallel distribution)
    }

    /**
     * @brief assigns a scalar to a vector
     * @details rX = A
     * @param rX The vector considered
     * @param A The scalar considered
     */
    static void Set(
        VectorType& rX,
        const DataType A
        )
    {
        rX.putScalar(static_cast<ST>(A));
    }

    /**
     * @brief Resizes a matrix
     * @param rA The matrix to be resized
     * @param m The new number of rows
     * @param n The new number of columns
     */
    static void Resize(
        MatrixType& rA,
        const SizeType m,
        const SizeType n
        )
    {
        KRATOS_ERROR << "Resize is not defined for Trilinos Sparse Matrix" << std::endl;
    }

    /**
     * @brief Resizes a vector
     * @param rX The vector to be resized
     * @param n The new size
     */
    static void Resize(
        VectorType& rX,
        const SizeType n
        )
    {
        KRATOS_ERROR << "Resize is not defined for a reference to Trilinos Vector - need to use the version passing a Pointer" << std::endl;
    }

    /**
     * @brief Resizes a vector
     * @param pA The pointer to the vector to be resized
     * @param n The new size
    */
    static void Resize(
        VectorPointerType pX,
        const SizeType n
        )
    {
        //KRATOS_ERROR_IF(pX != Teuchos::null) << "Trying to resize a null pointer" << std::endl;
        int global_elems = n;
        auto map = Teuchos::rcp(new MapType(0, 0, pX->getMap()->getComm()));
        VectorPointerType pNewEmptyX = Teuchos::rcp(new VectorType(map));
        pX.swap(pNewEmptyX);
    }

    /**
     * @brief Clears a matrix
     * @param pA The pointer to the matrix to be cleared
     */
    static void Clear(MatrixPointerType pA)
    {
        if(pA != Teuchos::null) {
            int global_elems = 0;
            auto map = Teuchos::rcp(new MapType(0, 0, pA->getMap()->getComm()));
            GraphPointerType graph = Teuchos::rcp(new GraphType(map, map, 0));
            MatrixPointerType pNewEmptyA = Teuchos::rcp(new MatrixType(graph));
            pA.swap(pNewEmptyA);
        }
    }

    /**
     * @brief Clears a vector
     * @param pX The pointer to the vector to be cleared
     */
    static void Clear(VectorPointerType pX)
    {
        if(pX != Teuchos::null) {
            int global_elems = 0;
            auto map = Teuchos::rcp(new MapType(0, 0, pX->getMap()->getComm()));
            VectorPointerType pNewEmptyX = Teuchos::rcp(new VectorType(map));
            pX.swap(pNewEmptyX);
        }
    }

    /**
     * @brief Sets a matrix to zero
     * @param rX The matrix to be set
     */
    inline static void SetToZero(MatrixType& rA)
    {
        // Set all values in the matrix to zero.
        rA.setAllToScalar(0.0);

        // Finalize the fill of the matrix if needed.
        rA.fillComplete();
    }

    /**
     * @brief Sets a vector to zero
     * @param rX The vector to be set
     */
    inline static void SetToZero(VectorType& rX)
    {
        rX.putScalar(static_cast<ST>(0));
    }

    /// TODO: creating the the calculating reaction version
    // 	template<class TOtherMatrixType, class TEquationIdVectorType>

    /**
     * @brief Assembles the LHS of the system
     * @param rA The LHS matrix
     * @param rLHSContribution The contribution to the LHS
     * @param rEquationId The equation ids
     */
    inline static void AssembleLHS(
        MatrixType& rA,
        const Matrix& rLHSContribution,
        const std::vector<std::size_t>& rEquationId
        )
    {
        const std::size_t system_size = rA.getGlobalNumRows();

        // Count active indices
        std::vector<LO> indices;
        for (std::size_t i = 0; i < rEquationId.size(); ++i) {
            if (rEquationId[i] < system_size) {
                indices.push_back(static_cast<LO>(rEquationId[i]));
            }
        }

        if (!indices.empty()) {
            // Fill Tpetra local arrays
            std::vector<ST> values(indices.size() * indices.size(), 0.0);

            for (std::size_t i = 0; i < indices.size(); ++i) {
                for (std::size_t j = 0; j < indices.size(); ++j) {
                    values[i * indices.size() + j] = rLHSContribution(i, j);
                }
            }

            // Insert or sum into global values
            Teuchos::ArrayView<LO> indices_view(indices);
            Teuchos::ArrayView<ST> values_view(values);

            const int ierr = rA.sumIntoGlobalValues(indices_view, indices_view, values_view);
            KRATOS_ERROR_IF(ierr != 0) << "Tpetra failure found" << std::endl;
        }
    }

    //***********************************************************************
    /// TODO: creating the the calculating reaction version
    // 	template<class TOtherVectorType, class TEquationIdVectorType>

    /**
     * @brief Assembles the RHS of the system
     * @param rb The RHS vector
     * @param rRHSContribution The RHS contribution
     * @param rEquationId The equation ids
     */
    inline static void AssembleRHS(
        VectorType& rb,
        const Vector& rRHSContribution,
        const std::vector<std::size_t>& rEquationId
        )
    {
        const std::size_t system_size = rb->getGlobalLength();

        // Count active indices
        std::vector<LO> indices;
        for (std::size_t i = 0; i < rEquationId.size(); ++i) {
            if (rEquationId[i] < system_size) {
                indices.push_back(static_cast<LO>(rEquationId[i]));
            }
        }

        if (!indices.empty()) {
            // Fill Tpetra values array
            std::vector<ST> values(indices.size(), 0.0);

            for (std::size_t i = 0; i < indices.size(); ++i) {
                values[i] = rRHSContribution[i];
            }

            // Insert or sum into global values
            Teuchos::ArrayView<LO> indices_view(indices);
            Teuchos::ArrayView<ST> values_view(values);

            const int ierr = rb->sumIntoGlobalValues(indices_view, values_view);
            KRATOS_ERROR_IF(ierr != 0) << "Tpetra failure found" << std::endl;
        }
    }

    /**
     * @brief This function returns if we are in a distributed system
     * @return True if we are in a distributed system, false otherwise (always true in this case)
     */
    inline static constexpr bool IsDistributed()
    {
        return true;
    }

    /**
     * @brief Returns a list of the fastest direct solvers.
     * @details This function returns a vector of strings representing the names of the fastest direct solvers. The order of the solvers in the list may need to be updated and reordered depending on the size of the equation system.
     * @return A vector of strings containing the names of the fastest direct solvers.
     */
    inline static std::vector<std::string> FastestDirectSolverList()
    {
        // May need to be updated and reordered. In fact I think it depends of the size of the equation system
        std::vector<std::string> faster_direct_solvers({
            "mumps2",         // Amesos2 (if compiled with MUMPS-support)
            "mumps",          // Amesos (if compiled with MUMPS-support)
            "super_lu_dist2", // Amesos2 SuperLUDist (if compiled with MPI-support)
            "super_lu_dist",  // Amesos SuperLUDist (if compiled with MPI-support)
            "amesos2",        // Amesos2
            "amesos",         // Amesos
            "klu2",           // Amesos2 KLU
            "klu",            // Amesos KLU
            "basker"          // Amesos2 Basker
        });
        return faster_direct_solvers;
    }

    /**
     * @brief This function returns a value from a given vector according to a given index
     * @param rX The vector from which values are to be gathered
     * @param I The index of the value to be gathered
     * @return The value of the vector corresponding to the index I
     */
    inline static double GetValue(
        const VectorType& rX,
        const std::size_t I
        )
    {
        // Get the local index corresponding to the global index `I`
        auto map = rX.getMap();
        IndexType localIndex = map->getLocalElement(static_cast<IndexType>(I));

        // Index must be local to this proc
        KRATOS_ERROR_IF(localIndex == Tpetra::Details::OrdinalTraits<IndexType>::invalid()) << " non-local id: " << I << "." << std::endl;

        // Get the value at the specified local index
        return rX.getLocalElement(localIndex);
    }

    /**
     * @brief This function gathers the values of a given vector according to a given index array
     * @param rX The vector from which values are to be gathered
     * @param IndexArray The array containing the indices of the values to be gathered
     * @param pValues The array containing the gathered values
     */
    static void GatherValues(
        const VectorType& rX,
        const std::vector<int>& IndexArray,
        double* pValues
        )
    {
        KRATOS_TRY

        // Get the total size of the index array
        const std::size_t tot_size = IndexArray.size();

        // Create a Map with the desired indices
        Teuchos::ArrayView<const IndexType> indexArrayView(IndexArray.data(), IndexArray.size());
        MapPointerType dof_update_map = Tpetra::createNonContigMapWithNode<IndexType, IndexType>(indexArrayView, rX.getMap()->getComm());

        // Define the Importer
        Tpetra::Import<IndexType, IndexType> importer(dof_update_map, rX.getMap());

        // Create a temporary vector to gather the values
        VectorType temp(dof_update_map);

        // Import the values from rX into the temp vector
        temp.doImport(rX, importer, Tpetra::INSERT);

        // Extract the values from the temp vector
        temp.get1dCopy(Teuchos::ArrayView<double>(pValues, tot_size));

        // Synchronize processes
        rX.getMap()->getComm()->barrier();

        KRATOS_CATCH("")
    }

    /**
     * @brief Read a matrix from a MatrixMarket file
     * @param rFileName The name of the file to read
     * @param rComm The MPI communicator
     * @return The matrix read from the file
     */
    MatrixPointerType ReadMatrixMarket(
        const std::string FileName,
        CommunicatorType& rComm
        )
    {
        KRATOS_TRY

    //     // Load the matrix using Tpetra's MatrixMarket reader
    //     RCP<const MapType> map;
    //     RCP<MatrixType> matrix = Tpetra::MatrixMarket::Reader<MatrixType>::readSparseFile(FileName, map, rComm);

    //     KRATOS_ERROR_IF(matrix.is_null()) << "Error thrown while reading Matrix Market file " << FileName;

    //     rComm.barrier();

    //     // Create a copy of the matrix
    //     MatrixPointerType paux = Teuchos::cp(new MatrixType(matrix->getRowMap(), matrix->getColMap(), matrix->getNodeNumEntries()));

    //     // Copy values from original matrix to the new matrix
    //     for (std::size_t i = 0; i < matrix->getNodeNumRows(); ++i) {
    //         const auto globalRow = matrix->getRowMap()->getGlobalElement(i);

    //         Teuchos::ArrayView<const ST> values;
    //         Teuchos::ArrayView<const LO> indices;

    //         matrix->getLocalRowView(i, indices, values);

    //         paux->replaceGlobalValues(globalRow, indices, values);
    //     }

    //     // Assemble the matrix
    //     paux->fillComplete();

    //     return paux;

        KRATOS_ERROR << "ReadMatrixMarket is not implemented for Trilinos TPetra" << std::endl;
        return CreateEmptyMatrixPointer();

        KRATOS_CATCH("");
    }

    /**
     * @brief Read a vector from a MatrixMarket file
     * @param rFileName The name of the file to read
     * @param pComm The MPI communicator
     * @param N The size of the vector
     */
    VectorPointerType ReadMatrixMarketVector(
        const std::string& rFileName,
        CommunicatorPointerType pComm,
        const int N
        )
    {
        KRATOS_TRY

    //     // Create a Tpetra::Map
    //     MapPointerType my_map = Tpetra::createUniformContigMapWithNode<GO, GO>(N, pComm);

    //     // Create an empty Tpetra vector
    //     VectorPointerType final_vector = Teuchos::rcp(new VectorType(my_map));

    //     // Read the vector from the Matrix Market file
    //     const int error_code = Tpetra::MatrixMarket::Reader<VectorType>::readVectorFile(rFileName, final_vector);
    //     KRATOS_ERROR_IF(error_code != 0) << "Error thrown while reading Matrix Market Vector file: " << rFileName << " error code is: " << error_code << std::endl;

    //     // Ensure all processes are synchronized
    //     pComm->barrier();

    //     // Global assemble the vector (if necessary)
    //     final_vector->doImport(*final_vector, Tpetra::Import<GO, GO>(my_map, my_map), Tpetra::INSERT);
    //     final_vector->globalAssemble();

    //     return final_vector;

        KRATOS_ERROR << "ReadMatrixMarketVector is not implemented for Trilinos TPetra" << std::endl;
        return CreateEmptyVectorPointer();

        KRATOS_CATCH("");
    }

    /**
    * @brief Generates a graph combining the graphs of two matrices
    * @param rA The first matrix
    * @param rB The second matrix
    */
    static GraphPointerType CombineMatricesGraphs(
        const MatrixType& rA,
        const MatrixType& rB
        )
    {
        // Row maps must be the same
        KRATOS_ERROR_IF(!rA.getRowMap()->isSameAs(*rB.getRowMap())) << "Row maps are not compatible" << std::endl;

        // Getting the graphs
        auto p_graph_a = rA.getCrsGraph();
        auto p_graph_b = rB.getCrsGraph();

        // Assuming local indexes
        KRATOS_ERROR_IF(!p_graph_a->isLocallyIndexed() || !p_graph_b->isLocallyIndexed()) << "Graphs indexes must be local" << std::endl;

        // Some definitions
        Teuchos::ArrayView<const int> cols_a, cols_b; // Column indices of row non-zero values (local indices should be 'int')
        const bool same_col_map = rA.getColMap()->isSameAs(*rB.getColMap());
        Teuchos::RCP<GraphType> graph;

        if (same_col_map) {
            // Same column map: create the graph with row and column maps from matrix A
            graph = Teuchos::rcp(new GraphType(rA.getRowMap(), rA.getColMap(), 1000));
        } else {
            // Different column maps: create the graph only with the row map from matrix A
            auto map = Teuchos::rcp(new MapType(0, 0, rA.getMap()->getComm()));
            graph = Teuchos::rcp(new GraphType(rA.getRowMap(), map, 1000));
        }

        const auto numLocalRows = p_graph_a->getRowMap()->getNodeNumElements();

        if (same_col_map) {
            // Combine graphs for the same column map using local indexing
            for (std::size_t i = 0; i < numLocalRows; ++i) {
                std::unordered_set<int> combined_indexes;

                // First graph
                p_graph_a->getLocalRowView(i, cols_a);
                combined_indexes.insert(cols_a.begin(), cols_a.end());

                // Second graph
                p_graph_b->getLocalRowView(i, cols_b);
                combined_indexes.insert(cols_b.begin(), cols_b.end());

                // Vector equivalent
                std::vector<int> combined_indexes_vector(combined_indexes.begin(), combined_indexes.end());

                // Adding to graph
                graph->insertLocalIndices(i, Teuchos::ArrayView<const int>(combined_indexes_vector));
            }
        } else {
            // Combine graphs with different column maps using global indexing
            for (std::size_t i = 0; i < numLocalRows; ++i) {
                const auto global_row_index = p_graph_a->getRowMap()->getGlobalElement(i);
                std::vector<long long int> combined_indexes_vector; // Global indices should be 'long long int'

                // First graph
                p_graph_a->getLocalRowView(i, cols_a);
                combined_indexes_vector.reserve(cols_a.size());
                for (long int j = 0; j < cols_a.size(); ++j) {
                    combined_indexes_vector.push_back(p_graph_a->getColMap()->getGlobalElement(cols_a[j]));
                }

                // Adding to graph
                graph->insertGlobalIndices(global_row_index, Teuchos::ArrayView<const long long int>(combined_indexes_vector));
            }

            for (std::size_t i = 0; i < numLocalRows; ++i) {
                const auto global_row_index = p_graph_b->getRowMap()->getGlobalElement(i);
                std::vector<long long int> combined_indexes_vector; // Global indices should be 'long long int'

                // Second graph
                p_graph_b->getLocalRowView(i, cols_b);
                combined_indexes_vector.reserve(cols_b.size());
                for (long int j = 0; j < cols_b.size(); ++j) {
                    combined_indexes_vector.push_back(p_graph_b->getColMap()->getGlobalElement(cols_b[j]));
                }

                // Adding to graph
                graph->insertGlobalIndices(global_row_index, Teuchos::ArrayView<const long long int>(combined_indexes_vector));
            }
        }

        // Finalizing graph construction
        graph->fillComplete();
        KRATOS_ERROR_IF(!graph->isFillComplete()) << "Tpetra graph fillComplete failed" << std::endl;

        return graph;
    }

    /**
     * @brief Copy values from one matrix to another
     * @details It is assumed that the sparsity of both matrices is compatible
     * @param rA The matrix where assigning values
     * @param rB The matrix to be copied
     */
    static void CopyMatrixValues(
        MatrixType& rA,
        const MatrixType& rB
        )
    {
        // Cleaning destination matrix
        SetToZero(rA);

        // Resume fill on the destination matrix
        rA.resumeFill();

        // Row maps must be the same
        const bool same_col_map = rA.getColMap()->isSameAs(*rB.getColMap());

        // Getting the graph of the source matrix
        auto p_graph_b = rB.getCrsGraph();

        // Copy values from rB to rA
        std::size_t num_entries; // Number of non-zero entries in rB matrix
        ColumnViewType cols;     // Column indices of row non-zero values (rB matrix)
        ValueViewType vals;      // Row non-zero values (rB matrix)

        if (same_col_map) {
            for (LO i = 0; i < static_cast<LO>(rB.getNodeNumRows()); ++i) {  // `i` should be a LocalOrdinal (LO)
                const auto global_row_index = rB.getRowMap()->getGlobalElement(i);
                rB.getGlobalRowCopy(global_row_index, cols, vals, num_entries);

                // Convert global column indices to local column indices
                Teuchos::Array<LO> local_cols(num_entries);
                for (size_t j = 0; j < num_entries; ++j) {
                    local_cols[j] = rA.getColMap()->getLocalElement(cols[j]);
                }

                // Sum values into local matrix using local row and column indices
                rA.sumIntoLocalValues(i, Teuchos::ArrayView<const LO>(local_cols), Teuchos::ArrayView<const ST>(vals.data(), num_entries));
            }
        } else {
            for (LO i = 0; i < static_cast<LO>(rB.getNodeNumRows()); ++i) {
                const auto global_row_index = rB.getRowMap()->getGlobalElement(i);
                rB.getGlobalRowCopy(global_row_index, cols, vals, num_entries);

                Teuchos::Array<GO> global_cols(num_entries);
                for (std::size_t j = 0; j < num_entries; ++j) {
                    global_cols[j] = p_graph_b->getColMap()->getGlobalElement(cols[j]);
                }

                // Sum values into global matrix using global row and column indices
                rA.sumIntoGlobalValues(global_row_index, Teuchos::ArrayView<const GO>(global_cols), Teuchos::ArrayView<const ST>(vals.data(), num_entries));
            }
        }

        // Complete fill of the matrix
        rA.fillComplete();
    }

    /**
     * @brief This method checks and corrects the zero diagonal values
     * @details This method returns the scale norm considering scaling the diagonal
     * @param rProcessInfo The problem process info
     * @param rA The LHS matrix
     * @param rb The RHS vector
     * @param ScalingDiagonal The type of scaling diagonal considered
     * @return The scale norm
     */
    static double CheckAndCorrectZeroDiagonalValues(
        const ProcessInfo& rProcessInfo,
        MatrixType& rA,
        VectorType& rb,
        const SCALING_DIAGONAL ScalingDiagonal = SCALING_DIAGONAL::NO_SCALING
        )
    {
        KRATOS_TRY

        // Define zero value tolerance
        const double zero_tolerance = std::numeric_limits<double>::epsilon();

        // The diagonal considered
        const double scale_factor = GetScaleNorm(rProcessInfo, rA, ScalingDiagonal);

        auto localMatrix = rA.getLocalMatrix();
        auto rowMap = rA.getRowMap();
        auto colMap = rA.getColMap();
        auto localRhs = rb.getLocalViewHost();

        for (int i = 0; i < localMatrix.numRows(); ++i) {
            auto localRow = localMatrix.row(i);
            const auto row_gid = rowMap->getGlobalElement(i);
            bool empty = true;
            int j;
            for (j = 0; j < localRow.length; ++j) {
                const auto col_gid = colMap->getGlobalElement(localRow.colidx(j));
                // Check diagonal value
                if (col_gid == row_gid) {
                    if (std::abs(localRow.value(j)) > zero_tolerance) {
                        empty = false;
                    }
                    break;
                }
            }

            // If diagonal empty assign scale factor
            if (empty) {
                const int row_gid_int = static_cast<int>(row_gid);  // Casting to int
                localMatrix.replaceValues(i, &row_gid_int, 1, &scale_factor, false, true);
                localRhs(i, 0) = 0.0;
            }
        }

        return scale_factor;

        KRATOS_CATCH("")
    }

    /**
     * @brief This method returns the scale norm considering for scaling the diagonal
     * @param rProcessInfo The problem process info
     * @param rA The LHS matrix
     * @param ScalingDiagonal The type of scaling diagonal considered
     * @return The scale norm
     */
    static double GetScaleNorm(
        const ProcessInfo& rProcessInfo,
        const MatrixType& rA,
        const SCALING_DIAGONAL ScalingDiagonal = SCALING_DIAGONAL::NO_SCALING
        )
    {
        KRATOS_TRY

        switch (ScalingDiagonal) {
            case SCALING_DIAGONAL::NO_SCALING:
                return 1.0;
            case SCALING_DIAGONAL::CONSIDER_PRESCRIBED_DIAGONAL: {
                KRATOS_ERROR_IF_NOT(rProcessInfo.Has(BUILD_SCALE_FACTOR)) << "Scale factor not defined at process info" << std::endl;
                return rProcessInfo.GetValue(BUILD_SCALE_FACTOR);
            }
            case SCALING_DIAGONAL::CONSIDER_NORM_DIAGONAL:
                return GetDiagonalNorm(rA)/static_cast<double>(Size1(rA));
            case SCALING_DIAGONAL::CONSIDER_MAX_DIAGONAL:
                return GetMaxDiagonal(rA);
            default:
                return GetMaxDiagonal(rA);
        }

        KRATOS_CATCH("");
    }

    /**
    * @brief This method returns the diagonal norm considering for scaling the diagonal
    * @param rA The LHS matrix
    * @return The diagonal norm
    */
    static double GetDiagonalNorm(const MatrixType& rA)
    {
        KRATOS_TRY

        // Create a vector to store the diagonal
        VectorType diag(rA.getRowMap());

        // Extract the diagonal entries
        rA.getLocalDiagCopy(diag);

        // Get the local view of the diagonal
        auto diagLocalView = diag.getLocalViewHost(Tpetra::Access::ReadOnly);

        // Compute the local sum of squares of the diagonal
        ST localSumOfSquares = Teuchos::ScalarTraits<ST>::zero();  // Initialize to 0

        auto numLocalEntries = diag.getLocalLength();
        ST value = 0.0;
        for (std::size_t i = 0; i < numLocalEntries; ++i) {
            value = diagLocalView(i, 0);
            localSumOfSquares += value * value;  // Sum of squares
        }

        // Perform a global reduction to sum the squares across all processes
        ST globalSumOfSquares = 0.0;
        Teuchos::reduceAll(*rA.getMap()->getComm(), Teuchos::REDUCE_SUM, localSumOfSquares, Teuchos::outArg(globalSumOfSquares));

        // Compute the two-norm by taking the square root of the global sum of squares
        return std::sqrt(globalSumOfSquares);

        KRATOS_CATCH("");
    }

    /**
    * @brief This method returns the diagonal max value
    * @param rA The LHS matrix
    * @return The diagonal max value
    */
    static double GetAveragevalueDiagonal(const MatrixType& rA)
    {
        KRATOS_TRY

        return 0.5 * (GetMaxDiagonal(rA) + GetMinDiagonal(rA));

        KRATOS_CATCH("");
    }

    /**
    * @brief This method returns the diagonal max value
    * @param rA The LHS matrix
    * @return The diagonal max value
    */
    static double GetMaxDiagonal(const MatrixType& rA)
    {
        KRATOS_TRY

        // Create a vector to store the diagonal
        VectorType diag(rA.getRowMap());

        // Extract the diagonal entries
        rA.getLocalDiagCopy(diag);

        // Get the local view of the diagonal
        auto diagLocalView = diag.getLocalViewHost(Tpetra::Access::ReadOnly);

        // Find the local maximum value
        ST localMax = Teuchos::ScalarTraits<ST>::zero();  // Initialize to 0

        auto numLocalEntries = diag.getLocalLength();
        ST value = 0.0;
        for (std::size_t i = 0; i < numLocalEntries; ++i) {
            value = diagLocalView(i, 0);
            if (value > localMax) {
                localMax = value;
            }
        }

        // Perform a global reduction to find the global maximum
        double globalMax = 0.0;
        Teuchos::reduceAll(*rA.getMap()->getComm(), Teuchos::REDUCE_MAX, localMax, Teuchos::outArg(globalMax));

        return globalMax;

        KRATOS_CATCH("");
    }

    /**
    * @brief This method returns the diagonal min value
    * @param rA The LHS matrix
    * @return The diagonal min value
    */
    static double GetMinDiagonal(const MatrixType& rA)
    {
        KRATOS_TRY

        // Create a vector to store the diagonal
        VectorType diag(rA.getRowMap());

        // Extract the diagonal entries
        rA.getLocalDiagCopy(diag);

        // Get the local view of the diagonal
        auto diagLocalView = diag.getLocalViewHost(Tpetra::Access::ReadOnly);

        // Find the local minimum value
        ST localMin = Teuchos::ScalarTraits<ST>::rmax();  // Initialize to the max possible value

        auto numLocalEntries = diag.getLocalLength();
        ST value = 0.0;
        for (std::size_t i = 0; i < numLocalEntries; ++i) {
            value = diagLocalView(i, 0);
            if (value < localMin) {
                localMin = value;
            }
        }

        // Perform a global reduction to find the global minimum
        double globalMin = 0.0;
        Teuchos::reduceAll(*rA.getMap()->getComm(), Teuchos::REDUCE_MIN, localMin, Teuchos::outArg(globalMin));

        return globalMin;

        KRATOS_CATCH("");
    }

   /**
    * @brief Check if the TrilinosSpaceExperimental is distributed.
    * @details This static member function checks whether the TrilinosSpaceExperimental is distributed or not.
    * @return True if the space is distributed, false otherwise.
    */
    static constexpr bool IsDistributedSpace()
    {
        return true;
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /**
     * @brief Turn back information as a string.
     * @return Info as a string.
     */
    virtual std::string Info() const
    {
        return "TrilinosSpaceExperimental";
    }

    /**
     * @brief Print information about this object.
     * @param rOStream The output stream to print on.
     */
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "TrilinosSpaceExperimental";
    }

    /**
     * @brief Print object's data.
     * @param rOStream The output stream to print on.
     */
    virtual void PrintData(std::ostream& rOStream) const
    {
    }

    // /**
    //  * @brief Writes a matrix to a file in MatrixMarket format
    //  * @param pFileName The name of the file to be written
    //  * @param rM The matrix to be written
    //  * @param Symmetric If the matrix is symmetric
    //  * @return True if the file was successfully written, false otherwise
    //  */
    // template< class TOtherMatrixType >
    // static bool WriteMatrixMarketMatrix(
    //     const char* pFileName,
    //     const TOtherMatrixType& rM,
    //     const bool Symmetric
    //     )
    // {
    //     // the argument "Symmetric" does not have an effect for Trilinos => needed for compatibility with other Spaces
    //     KRATOS_TRY;
    //     return EpetraExt::RowMatrixToMatrixMarketFile(pFileName, rM); // Returns 0 if no error, -1 if any problems with file system.
    //     KRATOS_CATCH("");
    // }

    // /**
    //  * @brief Writes a vector to a file in MatrixMarket format
    //  * @param pFileName The name of the file to be written
    //  * @param rV The vector to be written
    //  * @return True if the file was successfully written, false otherwise
    //  */
    // template< class VectorType >
    // static bool WriteMatrixMarketVector(
    //     const char* pFileName,
    //     const VectorType& rV
    //     )
    // {
    //     KRATOS_TRY;
    //     return EpetraExt::MultiVectorToMatrixMarketFile(pFileName, rV);
    //     KRATOS_CATCH("");
    // }

    // /**
    //  * @brief Creates a new dof updater
    //  * @return The new dof updater
    //  */
    // static DofUpdaterPointerType CreateDofUpdater()
    // {
    //     DofUpdaterType tmp;
    //     return tmp.Create();
    // }

    ///@}
private:
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    TrilinosSpaceExperimental & operator=(TrilinosSpaceExperimental const& rOther);

    /// Copy constructor.
    TrilinosSpaceExperimental(TrilinosSpaceExperimental const& rOther);

}; // Class TrilinosSpaceExperimental

///@}

} // namespace Kratos.