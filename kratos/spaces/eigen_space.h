//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes
#include <numeric>

// External includes

// Project includes
#include "includes/process_info.h"
#include "includes/default_interface.h"
#include "includes/matrix_market_interface.h"
#include "utilities/dof_updater.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

/**
 * @brief Diagonal scaling options of the builder-and-solvers
 * @details It is the same enumeration the uBLAS space declares. The scale factor is 1 for NO_SCALING, the diagonal
 * norm divided by the number of rows for CONSIDER_NORM_DIAGONAL, the largest absolute diagonal value for
 * CONSIDER_MAX_DIAGONAL and the BUILD_SCALE_FACTOR stored in the process info for CONSIDER_PRESCRIBED_DIAGONAL.
 */
enum class SCALING_DIAGONAL {NO_SCALING = 0, CONSIDER_NORM_DIAGONAL = 1, CONSIDER_MAX_DIAGONAL = 2, CONSIDER_PRESCRIBED_DIAGONAL = 3};

/**
 * @brief Forward declaration of the EigenSpace class template
 * @details It is declared here so that the TEigenSparseSpace and TEigenDenseSpace aliases can refer to it. The
 * class is defined below.
 * @tparam TDataType Data type used in the Eigen space, the scalar type of the matrix and vector entries.
 * @tparam TMatrixType Matrix type used in the Eigen space.
 * @tparam TVectorType Vector type used in the Eigen space.
 */
template <class TDataType, class TMatrixType, class TVectorType>
class EigenSpace;

/// EigenSpace alias for sparse matrices. It uses the row-major compressed (CSR) matrix EigenCompressedMatrix, with the
/// CSR index type TIndexType (KratosEigenIndexType by default), and the dense vector EigenVector, both of scalar type
/// TDataType.
template <class TDataType, class TIndexType = KratosEigenIndexType>
using TEigenSparseSpace =
    EigenSpace<TDataType, EigenCompressedMatrix<TDataType, TIndexType>, EigenVector<TDataType>>;

/// EigenSpace alias for dense matrices. It uses the dense matrix EigenMatrix and the dense vector EigenVector, both of
/// scalar type TDataType.
template <class TDataType>
using TEigenDenseSpace =
    EigenSpace<TDataType, EigenMatrix<TDataType>, EigenVector<TDataType>>;

///@}
///@name Kratos Classes
///@{

/**
 * @brief A class template for handling data types, matrices, and vectors in an Eigen-based space
 * @details Mirrors the static API of UblasSpace over the Eigen backend types (EigenCompressedMatrix, EigenMatrix,
 * EigenVector), so it can be used as TSparseSpace/TDenseSpace template argument by the strategies,
 * builder-and-solvers and linear solvers. It belongs to the KratosCore group and its author is Vicente Mataix
 * Ferrandiz.
 * @tparam TDataType Data type used in the Eigen space, the scalar type of the matrix and vector entries.
 * @tparam TMatrixType Matrix type used in the Eigen space.
 * @tparam TVectorType Vector type used in the Eigen space.
 */
template<class TDataType, class TMatrixType, class TVectorType>
class EigenSpace
{
public:
    ///@name Type Definitions
    ///@{

    /**
     * @brief Pointer definition of EigenSpace
     * @details Defines the Pointer and SharedPointer (Kratos::shared_ptr), WeakPointer (Kratos::weak_ptr) and
     * UniquePointer (Kratos::unique_ptr) aliases of the class.
     */
    KRATOS_CLASS_POINTER_DEFINITION(EigenSpace);

    /// The data type considered. Alias of the TDataType template argument, the scalar type of the matrix and vector
    /// entries.
    using DataType = TDataType;

    /// The matrix type considered. Alias of the TMatrixType template argument, the type of the matrices of the space.
    using MatrixType = TMatrixType;

    /// The vector type considered. Alias of the TVectorType template argument, the type of the vectors of the space.
    using VectorType = TVectorType;

    /// The index type considered. Unsigned integer type (std::size_t) used to index vector entries and matrix rows and
    /// columns.
    using IndexType = std::size_t;

    /// The size type considered. Unsigned integer type (std::size_t) used for the sizes of the vectors and matrices.
    using SizeType = std::size_t;

    /// The pointer to the matrix type. Shared pointer (Kratos::shared_ptr) to a matrix of the space.
    using MatrixPointerType = typename Kratos::shared_ptr<TMatrixType>;

    /// The pointer to the vector type. Shared pointer (Kratos::shared_ptr) to a vector of the space.
    using VectorPointerType = typename Kratos::shared_ptr<TVectorType>;

    /// The DoF updater type. The DofUpdater specialization for this space, the utility that updates the values of the
    /// DoF variables after solving the system.
    using DofUpdaterType = DofUpdater<EigenSpace<TDataType, TMatrixType, TVectorType>>;

    /// The pointer to the DoF updater type. Unique pointer to the DoF updater, as returned by CreateDofUpdater.
    using DofUpdaterPointerType = typename DofUpdaterType::UniquePointer;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     * @details Creates the space. The class has no data members, so the constructor does nothing.
     */
    EigenSpace()
    {
    }

    /**
     * @brief Destructor
     * @details Virtual destructor. The class has no data members, so it does nothing.
     */
    virtual ~EigenSpace()
    {
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Creates a pointer to an empty matrix
     * @details Returns a shared pointer to a newly allocated matrix with zero rows and zero columns.
     * @return Shared pointer to the new empty matrix, with zero rows and zero columns.
     */
    static MatrixPointerType CreateEmptyMatrixPointer()
    {
        return MatrixPointerType(new TMatrixType(0, 0));
    }

    /**
     * @brief Creates a pointer to an empty vector
     * @details Returns a shared pointer to a newly allocated vector with zero entries.
     * @return Shared pointer to the new empty vector, with zero entries.
     */
    static VectorPointerType CreateEmptyVectorPointer()
    {
        return VectorPointerType(new TVectorType(0));
    }

    /**
     * @brief This method returns true if the pointer is null
     * @details The pointer considered is pPointer, which can be of any type convertible to bool (e.g. a raw or
     * smart pointer). It returns true if the pointer is null.
     * @tparam TPointerType Type of the pointer, any type convertible to bool (e.g. a raw or smart pointer).
     * @param pPointer Pointer to be checked.
     * @return True if the pointer is null, false otherwise.
     */
    template<class TPointerType>
    inline static bool IsNull(const TPointerType& pPointer)
    {
        return !static_cast<bool>(pPointer);
    }

    /**
     * @brief Returns the size of a vector
     * @details Number of entries of the vector rV.
     * @param rV Vector whose size is requested.
     * @return Number of entries of rV.
     */
    static IndexType Size(VectorType const& rV)
    {
        return static_cast<IndexType>(rV.size());
    }

    /**
     * @brief Returns the number of rows of a matrix
     * @details Number of rows of the matrix rM.
     * @param rM Matrix whose number of rows is requested.
     * @return Number of rows of rM.
     */
    static IndexType Size1(MatrixType const& rM)
    {
        return rM.size1();
    }

    /**
     * @brief Returns the number of columns of a matrix
     * @details Number of columns of the matrix rM.
     * @param rM Matrix whose number of columns is requested.
     * @return Number of columns of rM.
     */
    static IndexType Size2(MatrixType const& rM)
    {
        return rM.size2();
    }

    /**
     * @brief Copies a column of a matrix into a vector, rX_i = rM_ij
     * @details It works on the dense Kratos::Matrix used for multi-column solves, as in UblasSpace. The column j of
     * rM is copied into rX, which is resized to the number of rows of rM if its size differs.
     * @tparam TColumnType Type of the vector that receives the column (it needs size, resize and the subscript
     * operator).
     * @param j Index of the column of rM to copy.
     * @param rM Dense matrix (Kratos::Matrix) whose column is copied.
     * @param rX Vector that receives the column; it is resized to the number of rows of rM if its size differs.
     */
    template<typename TColumnType>
    static void GetColumn(unsigned int j, Matrix& rM, TColumnType& rX)
    {
        if (rX.size() != rM.size1())
            rX.resize(rM.size1(), false);

        for (IndexType i = 0; i < rM.size1(); i++) {
            rX[i] = rM(i, j);
        }
    }

    /**
     * @brief Copies a vector into a column of a matrix, rM_ij = rX_i
     * @details It works on the dense Kratos::Matrix used for multi-column solves. The vector rX is written into the
     * column j of rM, so it must have at least as many entries as rM has rows. The matrix rM is not resized.
     * @tparam TColumnType Type of the vector copied into the column (it needs the subscript operator).
     * @param j Index of the column of rM to write.
     * @param rM Dense matrix (Kratos::Matrix) whose column is overwritten.
     * @param rX Vector copied into the column; it must have at least as many entries as rM has rows.
     */
    template<typename TColumnType>
    static void SetColumn(unsigned int j, Matrix& rM, TColumnType& rX)
    {
        for (IndexType i = 0; i < rM.size1(); i++) {
            rM(i, j) = rX[i];
        }
    }

    /**
     * @brief Copies a matrix, rY = rX
     * @details Assigns the matrix rX to the matrix rY.
     * @param rX Matrix to copy.
     * @param rY Matrix that receives the copy.
     */
    static void Copy(MatrixType const& rX, MatrixType& rY)
    {
        rY = rX;
    }

    /**
     * @brief Copies a vector, rY = rX
     * @details The vector rY is resized to the size of rX if needed. The entries are copied in parallel.
     * @param rX Vector to copy.
     * @param rY Vector that receives the copy; it is resized to the size of rX if needed.
     */
    static void Copy(VectorType const& rX, VectorType& rY)
    {
        const SizeType size = Size(rX);
        if (Size(rY) != size)
            rY.resize(size, false);
        const auto* x_data = rX.data().begin();
        auto* y_data = rY.data().begin();
        ParallelChunks(size, [=](const SizeType Begin, const SizeType End) {
            for (SizeType i = Begin; i < End; ++i)
                y_data[i] = x_data[i];
        });
    }

    /**
     * @brief Dot product of two vectors, rX * rY
     * @details The product is non-conjugated, as in UblasSpace. It is evaluated in parallel: the entries are split
     * in one contiguous chunk per thread and the partial dot products of the chunks are added up. The size of rX is
     * the number of entries considered.
     * @param rX First vector of the product; its size is the number of entries considered.
     * @param rY Second vector of the product, with at least as many entries as rX.
     * @return Non-conjugated dot product of rX and rY.
     */
    static TDataType Dot(VectorType const& rX, VectorType const& rY)
    {
        const SizeType size = Size(rX);
        const auto* x_data = rX.data().begin();
        const auto* y_data = rY.data().begin();
        const int n_chunks = ParallelUtilities::GetNumThreads();
        if (n_chunks == 1) {
            return PartialDot(x_data, y_data, SizeType(0), size);
        }
        return IndexPartition<int>(n_chunks, n_chunks).template for_each<SumReduction<TDataType>>([=](int Chunk) {
            const SizeType begin = (size * Chunk) / n_chunks;
            const SizeType end = (size * (Chunk + 1)) / n_chunks;
            return PartialDot(x_data, y_data, begin, end);
        });
    }

    /**
     * @brief Euclidean norm of a vector, ||rX||2
     * @details For real scalars it is the square root of the dot product of rX with itself. For complex scalars it
     * is the Euclidean norm sqrt(sum |x_i|^2), as the ublas norm_2.
     * @param rX Vector whose norm is computed.
     * @return Euclidean norm of rX.
     */
    static TDataType TwoNorm(VectorType const& rX)
    {
        if constexpr (Eigen::NumTraits<TDataType>::IsComplex) {
            return TDataType(rX.norm());
        } else {
            return std::sqrt(Dot(rX, rX));
        }
    }

    /**
     * @brief Frobenius norm of a dense matrix
     * @details Square root of the sum of the squares of all the entries of rA, evaluated in parallel by rows. For
     * complex scalars the squares are the plain powers std::pow(x, 2), without conjugation as in the uBLAS space, so
     * the result is then not the usual Frobenius norm.
     * @param rA Dense matrix whose norm is computed.
     * @return Square root of the sum of the squares of the entries of rA (the Frobenius norm for real scalars).
     */
    static TDataType TwoNorm(const EigenMatrix<TDataType>& rA) // Frobenius norm
    {
        const TDataType aux_sum = IndexPartition<IndexType>(rA.size1()).template for_each<SumReduction<TDataType>>([&](IndexType i) {
            TDataType row_sum = TDataType();
            for (IndexType j = 0; j < rA.size2(); ++j) {
                row_sum += std::pow(rA(i, j), 2);
            }
            return row_sum;
        });
        return std::sqrt(aux_sum);
    }

    /**
     * @brief TwoNorm of a (dense) Eigen expression, e.g. TwoNorm(y - x)
     * @details The exact template match avoids the ambiguity between the vector and dense-matrix overloads that a
     * converting expression argument would otherwise cause. The result is the Eigen norm of the expression: the
     * Euclidean norm for a vector and the Frobenius norm for a matrix.
     * @tparam TDerived Derived type of the Eigen dense expression.
     * @param rX Eigen dense expression (vector or matrix) whose norm is computed.
     * @return Eigen norm of rX: the Euclidean norm for a vector and the Frobenius norm for a matrix.
     */
    template<class TDerived>
    static TDataType TwoNorm(const Eigen::MatrixBase<TDerived>& rX)
    {
        return rX.norm();
    }

    /**
     * @brief Frobenius norm of a sparse matrix
     * @details Square root of the sum of the squares of the stored values of rA, evaluated in parallel by chunks of
     * the value array. For complex scalars the products are not conjugated (as in PartialDot), so the result is then
     * not the usual Frobenius norm. The matrix is first put in compressed mode if needed.
     * @tparam TIndexType Signed index type of the CSR storage of the sparse matrix.
     * @param rA Sparse matrix whose norm is computed.
     * @return Square root of the sum of the squares of the stored values of rA (the Frobenius norm for real scalars).
     */
    template<class TIndexType>
    static TDataType TwoNorm(const EigenCompressedMatrix<TDataType, TIndexType>& rA) // Frobenius norm
    {
        EnsureCompressed(rA);
        const auto* values = rA.valuePtr();
        const SizeType nnz = rA.nnz();
        const int n_chunks = ParallelUtilities::GetNumThreads();
        const TDataType aux_sum = IndexPartition<int>(n_chunks).template for_each<SumReduction<TDataType>>([=](int Chunk) {
            const SizeType begin = (nnz * Chunk) / n_chunks;
            const SizeType end = (nnz * (Chunk + 1)) / n_chunks;
            return PartialDot(values, values, begin, end);
        });
        return std::sqrt(aux_sum);
    }

    /**
     * @brief This method computes the Jacobi norm
     * @details Returns the Jacobi norm of the matrix rA, computed as the sum of the absolute values of its
     * off-diagonal entries. The sum is evaluated in parallel by rows.
     * @param rA Dense matrix whose Jacobi norm is computed.
     * @return Sum of the absolute values of the off-diagonal entries of rA.
     */
    static TDataType JacobiNorm(const EigenMatrix<TDataType>& rA)
    {
        return IndexPartition<IndexType>(rA.size1()).template for_each<SumReduction<TDataType>>([&](IndexType i) {
            TDataType row_sum = TDataType();
            for (IndexType j = 0; j < rA.size2(); ++j) {
                if (i != j) {
                    row_sum += std::abs(rA(i, j));
                }
            }
            return row_sum;
        });
    }

    /**
     * @brief This method computes the Jacobi norm of a sparse matrix
     * @details Returns the Jacobi norm of the matrix rA, computed as the sum of the absolute values of its stored
     * off-diagonal entries. The matrix is first put in compressed mode if needed, then the sum is evaluated in
     * parallel by rows over the CSR arrays.
     * @tparam TIndexType Signed index type of the CSR storage of the sparse matrix.
     * @param rA Sparse matrix whose Jacobi norm is computed.
     * @return Sum of the absolute values of the stored off-diagonal entries of rA.
     */
    template<class TIndexType>
    static TDataType JacobiNorm(const EigenCompressedMatrix<TDataType, TIndexType>& rA)
    {
        EnsureCompressed(rA);
        const auto* row_ptr = rA.outerIndexPtr();
        const auto* col_idx = rA.innerIndexPtr();
        const auto* values = rA.valuePtr();
        return IndexPartition<IndexType>(rA.size1()).template for_each<SumReduction<TDataType>>([&](IndexType i) {
            TDataType row_sum = TDataType();
            for (auto k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
                if (static_cast<IndexType>(col_idx[k]) != i) {
                    row_sum += std::abs(values[k]);
                }
            }
            return row_sum;
        });
    }

    /**
     * @brief Dense matrix-vector product, rY = rA * rX
     * @details The vector rY is resized to the number of rows of rA if its size differs. The product is assigned
     * assuming that rY does not alias rA or rX, so no temporary is created.
     * @param rA Dense matrix of the product.
     * @param rX Vector multiplied by rA.
     * @param rY Vector that receives the product; it is resized to the number of rows of rA if its size differs.
     */
    static void Mult(const EigenMatrix<TDataType>& rA, const VectorType& rX, VectorType& rY)
    {
        if (Size(rY) != rA.size1())
            rY.resize(rA.size1(), false);
        rY.noalias() = rA * rX;
    }

    /**
     * @brief Sparse matrix-vector product (SpMV), rY = rA * rX
     * @details Deterministic manual CSR loop (as in KratosSpace) rather than Eigen's internally parallelized
     * product, to avoid nested parallelism under the Kratos shared-memory parallelization and keep UblasSpace
     * behavior. The per-chunk work is delegated to a plain function taking the raw pointers as arguments (as
     * UblasSpace does): keeping the hot loop out of the OpenMP-outlined lambda lets the compiler optimize it fully,
     * which benchmarks at parity with the uBLAS implementation. The vector rY is resized to the number of rows of
     * rA if needed.
     * @tparam TIndexType Signed index type of the CSR storage of the sparse matrix.
     * @param rA Sparse matrix of the product.
     * @param rX Vector multiplied by rA.
     * @param rY Vector that receives the product; it is resized to the number of rows of rA if its size differs.
     */
    template<class TIndexType>
    static void Mult(const EigenCompressedMatrix<TDataType, TIndexType>& rA, const VectorType& rX, VectorType& rY)
    {
        const SizeType n_rows = rA.size1();
        if (Size(rY) != n_rows)
            rY.resize(n_rows, false);

        EnsureCompressed(rA);
        const auto* row_ptr = rA.outerIndexPtr();
        const auto* col_idx = rA.innerIndexPtr();
        const auto* values = rA.valuePtr();
        const auto* x_data = rX.data().begin();
        auto* y_data = rY.data().begin();
        const int n_chunks = ParallelUtilities::GetNumThreads();
        if (n_chunks == 1) {
            PartialSpMV(row_ptr, col_idx, values, x_data, y_data, SizeType(0), n_rows);
            return;
        }
        IndexPartition<int>(n_chunks, n_chunks).for_each([=](int Chunk) {
            const SizeType begin = (n_rows * Chunk) / n_chunks;
            const SizeType end = (n_rows * (Chunk + 1)) / n_chunks;
            PartialSpMV(row_ptr, col_idx, values, x_data, y_data, begin, end);
        });
    }

    // rY = rA^T * rX
    /**
     * @brief Transpose matrix-vector product, rY = rA^T * rX
     * @details The vector rY is resized to the number of columns of rA if its size differs. The product is assigned
     * assuming that rY does not alias rA or rX, so no temporary is created.
     * @tparam TOtherMatrixType Type of the matrix (it needs size2 and transpose, e.g. an Eigen matrix).
     * @param rA Matrix whose transpose multiplies rX.
     * @param rX Vector multiplied by the transpose of rA.
     * @param rY Vector that receives the product; it is resized to the number of columns of rA if its size differs.
     */
    template<class TOtherMatrixType>
    static void TransposeMult(const TOtherMatrixType& rA, const VectorType& rX, VectorType& rY)
    {
        if (Size(rY) != rA.size2())
            rY.resize(rA.size2(), false);
        rY.noalias() = rA.transpose() * rX;
    }

    /**
     * @brief Sparse transpose matrix-vector product, rY = rA^T * rX
     * @details A deterministic serial CSR scatter (y[col] += v * x[row]), exactly the access pattern of the uBLAS
     * axpy_prod(x, A, y) implementation (which is serial too). Eigen's expression path for this case does an extra
     * full setZero() pass over rY and re-tests isCompressed() per row inside its InnerIterator, which benchmarks
     * measurably slower. The vector rY is resized to the number of columns of rA if needed and zeroed first.
     * @tparam TIndexType Signed index type of the CSR storage of the sparse matrix.
     * @param rA Sparse matrix whose transpose multiplies rX.
     * @param rX Vector multiplied by the transpose of rA.
     * @param rY Vector that receives the product; it is resized to the number of columns of rA if needed and zeroed
     * first.
     */
    template<class TIndexType>
    static void TransposeMult(const EigenCompressedMatrix<TDataType, TIndexType>& rA, const VectorType& rX, VectorType& rY)
    {
        const SizeType n_rows = rA.size1();
        const SizeType n_cols = rA.size2();
        if (Size(rY) != n_cols)
            rY.resize(n_cols, false);

        EnsureCompressed(rA);
        const auto* row_ptr = rA.outerIndexPtr();
        const auto* col_idx = rA.innerIndexPtr();
        const auto* values = rA.valuePtr();
        const auto* x_data = rX.data().begin();
        auto* y_data = rY.data().begin();

        SetToZero(rY);
        for (SizeType i = 0; i < n_rows; ++i) {
            const TDataType x_i = x_data[i];
            const auto row_end = row_ptr[i + 1];
            for (auto k = row_ptr[i]; k < row_end; ++k) {
                y_data[col_idx[k]] += values[k] * x_i;
            }
        }
    }

    /**
     * @brief Degree of a row in the graph of a matrix
     * @details Number of stored entries of the row i of the matrix rA, taken from the CSR row pointers. The matrix
     * is first put in compressed mode if needed.
     * @param i Index of the row.
     * @param rA CSR matrix whose graph is considered; it is put in compressed mode if needed.
     * @return Number of stored entries of the row i of rA.
     */
    static inline SizeType GraphDegree(IndexType i, TMatrixType& rA)
    {
        EnsureCompressed(rA);
        const auto* row_ptr = rA.outerIndexPtr();
        return static_cast<SizeType>(row_ptr[i + 1] - row_ptr[i]);
    }

    /**
     * @brief Neighbors of a row in the graph of a matrix
     * @details The vector rNeighbors is cleared and filled with the column indices of the stored entries of the row
     * i of the matrix rA, which is first put in compressed mode if needed.
     * @param i Index of the row.
     * @param rA CSR matrix whose graph is considered; it is put in compressed mode if needed.
     * @param rNeighbors Vector that is cleared and filled with the column indices of the stored entries of the row i.
     */
    static inline void GraphNeighbors(IndexType i, TMatrixType& rA, std::vector<IndexType>& rNeighbors)
    {
        EnsureCompressed(rA);
        rNeighbors.clear();
        const auto* row_ptr = rA.outerIndexPtr();
        const auto* col_idx = rA.innerIndexPtr();
        for (auto k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
            rNeighbors.push_back(static_cast<IndexType>(col_idx[k]));
        }
    }

    //********************************************************************
    // checks if a multiplication is needed and tries to do otherwise
    /**
     * @brief Multiplies a vector by a scalar in place, rX *= A
     * @details Nothing is done if A is 1, as no multiplication is needed. Otherwise the entries are scaled in
     * parallel by chunks.
     * @param rX Vector scaled in place.
     * @param A Scalar factor.
     */
    static void InplaceMult(VectorType& rX, const double A)
    {
        if (A == 1.00)
            return;
        auto* x_data = rX.data().begin();
        ParallelChunks(Size(rX), [=](const SizeType Begin, const SizeType End) {
            for (SizeType i = Begin; i < End; ++i)
                x_data[i] *= A;
        });
    }

    //********************************************************************
    // ATTENTION it is assumed no aliasing between rX and rY
    // X = A*y;
    /**
     * @brief Assigns a scaled vector, rX = A * rY
     * @details The vector rX is resized to the size of rY if their sizes differ. It is assumed that there is no
     * aliasing between rX and rY. The cases A equal to 1 and A equal to -1 have a fast path (plain and negated
     * copy), as in UblasSpace. The entries are assigned in parallel by chunks.
     * @param rX Vector that receives the result; it is resized to the size of rY if their sizes differ.
     * @param A Scalar factor that multiplies rY.
     * @param rY Vector to scale and assign to rX.
     */
    static void Assign(VectorType& rX, const double A, const VectorType& rY)
    {
        const SizeType size = Size(rY);
        if (Size(rX) != size)
            rX.resize(size, false);

        const auto* y_data = rY.data().begin();
        auto* x_data = rX.data().begin();
        // A == +-1 fast paths, as in UblasSpace
        if (A == 1.00) {
            ParallelChunks(size, [=](const SizeType Begin, const SizeType End) {
                for (SizeType i = Begin; i < End; ++i)
                    x_data[i] = y_data[i];
            });
        } else if (A == -1.00) {
            ParallelChunks(size, [=](const SizeType Begin, const SizeType End) {
                for (SizeType i = Begin; i < End; ++i)
                    x_data[i] = -y_data[i];
            });
        } else {
            ParallelChunks(size, [=](const SizeType Begin, const SizeType End) {
                for (SizeType i = Begin; i < End; ++i)
                    x_data[i] = A * y_data[i];
            });
        }
    }

    //********************************************************************
    // ATTENTION it is assumed no aliasing between rX and rY
    // X += A*y;
    /**
     * @brief Adds a scaled vector, rX += A * rY
     * @details The vector rX is resized to the size of rY (without preserving its values) if their sizes differ. It
     * is assumed that there is no aliasing between rX and rY. The cases A equal to 1 and A equal to -1 have a fast
     * path (plain add and subtract), as in UblasSpace. The entries are updated in parallel by chunks.
     * @param rX Vector to which the scaled vector is added; it is resized to the size of rY if their sizes differ.
     * @param A Scalar factor that multiplies rY.
     * @param rY Vector to scale and add to rX.
     */
    static void UnaliasedAdd(VectorType& rX, const double A, const VectorType& rY)
    {
        const SizeType size = Size(rY);
        if (Size(rX) != size)
            rX.resize(size, false);

        const auto* y_data = rY.data().begin();
        auto* x_data = rX.data().begin();
        // A == +-1 fast paths, as in UblasSpace
        if (A == 1.00) {
            ParallelChunks(size, [=](const SizeType Begin, const SizeType End) {
                for (SizeType i = Begin; i < End; ++i)
                    x_data[i] += y_data[i];
            });
        } else if (A == -1.00) {
            ParallelChunks(size, [=](const SizeType Begin, const SizeType End) {
                for (SizeType i = Begin; i < End; ++i)
                    x_data[i] -= y_data[i];
            });
        } else {
            ParallelChunks(size, [=](const SizeType Begin, const SizeType End) {
                for (SizeType i = Begin; i < End; ++i)
                    x_data[i] += A * y_data[i];
            });
        }
    }

    //********************************************************************

    // Fused single-pass implementations (one parallel region, rZ/rY written
    // once instead of the Assign + UnaliasedAdd two-pass composition, which
    // touches the target three times). The generic multiply subsumes the
    // +-1 special cases at no cost here: the fused loop is one FMA either way.

    /**
     * @brief Fused scale and add into a third vector, rZ = (A * rX) + (B * rY)
     * @details Single-pass implementation: one parallel region with rZ written once, instead of the Assign +
     * UnaliasedAdd two-pass composition, which touches the target three times. The generic multiply subsumes the
     * +-1 special cases at no cost here, as the fused loop is one FMA either way. The vector rZ is resized to the
     * size of rX if needed.
     * @param A Scalar factor that multiplies rX.
     * @param rX First vector.
     * @param B Scalar factor that multiplies rY.
     * @param rY Second vector, with at least as many entries as rX.
     * @param rZ Vector that receives the result; it is resized to the size of rX if needed.
     */
    static void ScaleAndAdd(const double A, const VectorType& rX, const double B, const VectorType& rY, VectorType& rZ) // rZ = (A * rX) + (B * rY)
    {
        const SizeType size = Size(rX);
        if (Size(rZ) != size)
            rZ.resize(size, false);

        const auto* x_data = rX.data().begin();
        const auto* y_data = rY.data().begin();
        auto* z_data = rZ.data().begin();
        ParallelChunks(size, [=](const SizeType Begin, const SizeType End) {
            for (SizeType i = Begin; i < End; ++i)
                z_data[i] = A * x_data[i] + B * y_data[i];
        });
    }

    /**
     * @brief Fused scale and add in place, rY = (A * rX) + (B * rY)
     * @details Single-pass implementation: one parallel region with rY written once, instead of the Assign +
     * UnaliasedAdd two-pass composition. The vector rY is not resized and must have at least the size of rX.
     * @param A Scalar factor that multiplies rX.
     * @param rX Vector scaled and added.
     * @param B Scalar factor that multiplies the previous value of rY.
     * @param rY Vector updated in place; it must have at least the size of rX.
     */
    static void ScaleAndAdd(const double A, const VectorType& rX, const double B, VectorType& rY) // rY = (A * rX) + (B * rY)
    {
        const SizeType size = Size(rX);
        const auto* x_data = rX.data().begin();
        auto* y_data = rY.data().begin();
        ParallelChunks(size, [=](const SizeType Begin, const SizeType End) {
            for (SizeType i = Begin; i < End; ++i)
                y_data[i] = A * x_data[i] + B * y_data[i];
        });
    }

    /**
     * @brief Dot product of a row of a matrix and a vector, rA[i] * rX
     * @details The product is non-conjugated, as in UblasSpace. It accumulates rA(i, j) * rX[j] over all the
     * columns j of rA.
     * @param i Index of the row of rA.
     * @param rA Matrix whose row i is used.
     * @param rX Vector multiplied by the row.
     * @return Non-conjugated dot product of the row i of rA and rX.
     */
    static double RowDot(unsigned int i, MatrixType& rA, VectorType& rX)
    {
        TDataType accumulator = TDataType();
        for (IndexType j = 0; j < rA.size2(); ++j) {
            accumulator += rA(i, j) * rX[j];
        }
        return accumulator;
    }

    /**
     * @brief Sets an entry of a vector
     * @details Assigns the given value to the entry i of the vector rX.
     * @param rX Vector to modify.
     * @param i Index of the entry to set.
     * @param value Value assigned to the entry.
     */
    static void SetValue(VectorType& rX, IndexType i, TDataType value)
    {
        rX[i] = value;
    }

    /**
     * @brief Sets all the entries of a vector to a value, rX = A
     * @details Every entry of the vector rX is set to the value A, keeping its size.
     * @param rX Vector to modify.
     * @param A Value assigned to every entry.
     */
    static void Set(VectorType& rX, TDataType A)
    {
        rX.setConstant(A);
    }

    /**
     * @brief Resizes a matrix
     * @details The matrix rA is resized to m rows and n columns without preserving its previous content (the values
     * and, for a sparse matrix, the structure are discarded).
     * @param rA Matrix to resize.
     * @param m New number of rows.
     * @param n New number of columns.
     */
    static void Resize(MatrixType& rA, SizeType m, SizeType n)
    {
        rA.resize(m, n, false);
    }

    /**
     * @brief Resizes the matrix pointed to
     * @details The matrix pointed to by pA is resized to m rows and n columns without preserving its previous
     * content (the values and, for a sparse matrix, the structure are discarded).
     * @param pA Pointer to the matrix to resize.
     * @param m New number of rows.
     * @param n New number of columns.
     */
    static void Resize(MatrixPointerType& pA, SizeType m, SizeType n)
    {
        pA->resize(m, n, false);
    }

    /**
     * @brief Resizes a vector
     * @details The vector rX is resized to n entries without preserving its previous values.
     * @param rX Vector to resize.
     * @param n New number of entries.
     */
    static void Resize(VectorType& rX, SizeType n)
    {
        rX.resize(n, false);
    }

    /**
     * @brief Resizes the vector pointed to
     * @details The vector pointed to by pX is resized to n entries without preserving its previous values.
     * @param pX Pointer to the vector to resize.
     * @param n New number of entries.
     */
    static void Resize(VectorPointerType& pX, SizeType n)
    {
        pX->resize(n, false);
    }

    /**
     * @brief Clears the matrix pointed to
     * @details The matrix pointed to by pA is resized to zero rows and zero columns, discarding its content. The
     * pointer itself is left untouched.
     * @param pA Pointer to the matrix to clear.
     */
    static void Clear(MatrixPointerType& pA)
    {
        pA->resize(0, 0, false);
    }

    /**
     * @brief Clears the vector pointed to
     * @details The vector pointed to by pX is resized to zero entries, discarding its content. The pointer itself
     * is left untouched.
     * @param pX Pointer to the vector to clear.
     */
    static void Clear(VectorPointerType& pX)
    {
        pX->resize(0, false);
    }

    /**
     * @brief Resizes the data storage of a matrix
     * @details For a sparse matrix (one providing resizeNonZeros) the value storage is resized to m non-zero
     * entries keeping the graph arrays, and all the values are set to zero. For any other matrix type an error is
     * raised, as it is not implemented.
     * @tparam TOtherMatrixType Type of the matrix; only types providing resizeNonZeros (sparse matrices) are supported.
     * @param rA Matrix whose data storage is resized.
     * @param m New number of non-zero entries of the value storage.
     */
    template<class TOtherMatrixType>
    inline static void ResizeData(TOtherMatrixType& rA, SizeType m)
    {
        if constexpr (requires { rA.resizeNonZeros(m); }) {
            // Sparse: resize the value storage keeping the graph arrays
            rA.resizeNonZeros(m);
            auto values = rA.value_data();
            std::fill(values.begin(), values.end(), TDataType());
        } else {
            KRATOS_ERROR << "ResizeData is not implemented for this matrix type" << std::endl;
        }
    }

    /**
     * @brief Resizes the data storage of a vector
     * @details The vector rX is resized to m entries (the previous values are not preserved) and all of them are
     * set to zero.
     * @param rX Vector whose data storage is resized.
     * @param m New number of entries.
     */
    inline static void ResizeData(VectorType& rX, SizeType m)
    {
        rX.resize(m, false);
        rX.setZero();
    }

    /**
     * @brief Sets a matrix to zero
     * @details For a sparse matrix the stored values are set to zero in parallel but the sparsity pattern is kept,
     * as setZero() would drop the graph; the matrix is first put in compressed mode if needed. For a dense matrix
     * all the entries are set to zero.
     * @param rA Matrix to set to zero.
     */
    inline static void SetToZero(MatrixType& rA)
    {
        if constexpr (requires { rA.valuePtr(); }) {
            // Sparse: zero the values but keep the sparsity pattern
            // (setZero() would drop the graph)
            EnsureCompressed(rA);
            auto* values = rA.valuePtr();
            ParallelChunks(rA.nnz(), [=](const SizeType Begin, const SizeType End) {
                for (SizeType i = Begin; i < End; ++i)
                    values[i] = TDataType();
            });
        } else {
            rA.setZero();
        }
    }

    /**
     * @brief Sets a vector to zero
     * @details All the entries of the vector rX are set to zero in parallel by chunks.
     * @param rX Vector to set to zero.
     */
    inline static void SetToZero(VectorType& rX)
    {
        auto* x_data = rX.data().begin();
        ParallelChunks(Size(rX), [=](const SizeType Begin, const SizeType End) {
            for (SizeType i = Begin; i < End; ++i)
                x_data[i] = TDataType();
        });
    }

    /**
     * @brief Assembles a local LHS contribution into the global matrix
     * @details Adds the entries of the local matrix rLHSContribution to rA, using the vector rEquationId to map the
     * local rows and columns to the global ones. The entries whose row or column equation id is not lower than the
     * size of the system are skipped.
     * @tparam TOtherMatrixType Type of the local matrix (it needs size1 and the call operator).
     * @tparam TEquationIdVectorType Type of the equation id vector (it needs the subscript operator).
     * @param rA Global matrix to which the local contribution is added.
     * @param rLHSContribution Local matrix with the contribution to assemble.
     * @param rEquationId Equation ids that map the local rows and columns to the global ones.
     */
    template<class TOtherMatrixType, class TEquationIdVectorType>
    inline static void AssembleLHS(
        MatrixType& rA,
        TOtherMatrixType& rLHSContribution,
        TEquationIdVectorType& rEquationId
        )
    {
        const unsigned int system_size = rA.size1();
        const unsigned int local_size = rLHSContribution.size1();

        for (unsigned int i_local = 0; i_local < local_size; i_local++) {
            const unsigned int i_global = rEquationId[i_local];
            if (i_global < system_size) {
                for (unsigned int j_local = 0; j_local < local_size; j_local++) {
                    const unsigned int j_global = rEquationId[j_local];
                    if (j_global < system_size)
                        rA(i_global, j_global) += rLHSContribution(i_local, j_local);
                }
            }
        }
    }

    /**
     * @brief This method checks and corrects the zero diagonal values
     * @details This method returns the scale norm considering for scaling the diagonal, as computed by GetScaleNorm
     * from the problem process info rProcessInfo, the LHS matrix rA and the type of scaling diagonal considered,
     * ScalingDiagonal. Each row of rA without any stored entry larger than the machine epsilon (in absolute value)
     * gets its diagonal entry set to the scale norm, and the corresponding entry of the RHS vector rb is set to
     * zero. When the diagonal term of such a row is not defined in the matrix, a warning is printed, and an error is
     * also raised in debug builds. The matrix is first put in compressed mode if needed.
     * @param rProcessInfo Process info of the problem, it provides BUILD_SCALE_FACTOR for the prescribed diagonal.
     * @param rA LHS matrix to correct.
     * @param rb RHS vector, whose entries of the corrected rows are set to zero.
     * @param ScalingDiagonal Type of scaling diagonal considered (NO_SCALING by default).
     * @return Scale norm obtained from GetScaleNorm, the value set on the corrected diagonal entries.
     */
    static double CheckAndCorrectZeroDiagonalValues(
        const ProcessInfo& rProcessInfo,
        MatrixType& rA,
        VectorType& rb,
        const SCALING_DIAGONAL ScalingDiagonal = SCALING_DIAGONAL::NO_SCALING
        )
    {
        // The system size
        const std::size_t system_size = rA.size1();

        // The matrix data
        EnsureCompressed(rA);
        auto values = rA.value_data();
        const auto row_indices = rA.index1_data();
        const auto col_indices = rA.index2_data();

        // Define the iterators
        const auto it_col_indices_begin = col_indices.begin();

        // Define zero value tolerance
        const double zero_tolerance = std::numeric_limits<double>::epsilon();

        // The diagonal considered
        const double scale_factor = GetScaleNorm(rProcessInfo, rA, ScalingDiagonal);

        // Detect if there is a line of all zeros and set the diagonal to a 1 if this happens
        IndexPartition<std::size_t>(system_size).for_each([&](std::size_t Index) {
            bool empty = true;

            const std::size_t col_begin = row_indices[Index];
            const std::size_t col_end = row_indices[Index + 1];

            for (std::size_t j = col_begin; j < col_end; ++j) {
                if (std::abs(values[j]) > zero_tolerance) {
                    empty = false;
                    break;
                }
            }

            if (empty) {
                const auto it_col_indices_row_begin = it_col_indices_begin + col_begin;
                const auto it_col_indices_row_end = it_col_indices_begin + col_end;

                const auto lower = std::lower_bound(it_col_indices_row_begin, it_col_indices_row_end, static_cast<typename MatrixType::StorageIndex>(Index));
                const auto upper = std::upper_bound(it_col_indices_row_begin, it_col_indices_row_end, static_cast<typename MatrixType::StorageIndex>(Index));

                if (lower != upper) { // Index was found
                    values[std::distance(it_col_indices_begin, lower)] = scale_factor;
                } else {
                    KRATOS_DEBUG_ERROR << "Diagonal term (" << Index << ", " << Index << ") is not defined in the system matrix" << std::endl;
                    KRATOS_WARNING("EigenSpace") << "Diagonal term (" << Index << ", " << Index << ") is not defined in the system matrix" << std::endl;
                }
                rb[Index] = 0.0;
            }
        });

        return scale_factor;
    }

    /**
     * @brief This method returns the scale norm considering for scaling the diagonal
     * @details The scale norm is selected by the type of scaling diagonal considered, ScalingDiagonal: 1 for
     * NO_SCALING, the BUILD_SCALE_FACTOR of the problem process info rProcessInfo for CONSIDER_PRESCRIBED_DIAGONAL
     * (an error is raised if it is not defined), the diagonal norm of the LHS matrix rA divided by its number of
     * rows for CONSIDER_NORM_DIAGONAL, and the maximum diagonal value of rA for CONSIDER_MAX_DIAGONAL (also for any
     * other value). The scale norm is returned.
     * @param rProcessInfo Process info of the problem, it provides BUILD_SCALE_FACTOR for the prescribed diagonal.
     * @param rA LHS matrix of the system.
     * @param ScalingDiagonal Type of scaling diagonal considered (NO_SCALING by default).
     * @return Scale norm selected by ScalingDiagonal.
     */
    static double GetScaleNorm(
        const ProcessInfo& rProcessInfo,
        const MatrixType& rA,
        const SCALING_DIAGONAL ScalingDiagonal = SCALING_DIAGONAL::NO_SCALING
        )
    {
        switch (ScalingDiagonal) {
            case SCALING_DIAGONAL::NO_SCALING:
                return 1.0;
            case SCALING_DIAGONAL::CONSIDER_PRESCRIBED_DIAGONAL: {
                KRATOS_ERROR_IF_NOT(rProcessInfo.Has(BUILD_SCALE_FACTOR)) << "Scale factor not defined at process info" << std::endl;
                return rProcessInfo.GetValue(BUILD_SCALE_FACTOR);
            }
            case SCALING_DIAGONAL::CONSIDER_NORM_DIAGONAL:
                return GetDiagonalNorm(rA)/static_cast<double>(rA.size1());
            case SCALING_DIAGONAL::CONSIDER_MAX_DIAGONAL:
                return GetMaxDiagonal(rA);
            default:
                return GetMaxDiagonal(rA);
        }
    }

    /**
     * @brief This method returns the diagonal norm considering for scaling the diagonal
     * @details Returns the diagonal norm of the LHS matrix rA, i.e. the square root of the sum of the squared
     * absolute values of its diagonal entries (a row without a stored diagonal entry contributes zero). The matrix
     * is first put in compressed mode if needed, then the sum is evaluated in parallel by rows.
     * @param rA LHS matrix of the system.
     * @return Square root of the sum of the squared absolute values of the diagonal entries of rA.
     */
    static double GetDiagonalNorm(const MatrixType& rA)
    {
        EnsureCompressed(rA);
        const auto values = rA.value_data();
        const auto row_indices = rA.index1_data();
        const auto col_indices = rA.index2_data();

        const double diagonal_norm = IndexPartition<std::size_t>(Size1(rA)).template for_each<SumReduction<double>>([&](std::size_t Index) {
            const std::size_t col_begin = row_indices[Index];
            const std::size_t col_end = row_indices[Index + 1];
            for (std::size_t j = col_begin; j < col_end; ++j) {
                if (static_cast<std::size_t>(col_indices[j]) == Index) {
                    return std::pow(std::abs(values[j]), 2);
                }
            }
            return 0.0;
        });

        return std::sqrt(diagonal_norm);
    }

    /**
     * @brief This method returns the diagonal average value
     * @details Returns the diagonal average value of the LHS matrix rA, computed as the mean of its maximum and
     * minimum diagonal values (see GetMaxDiagonal and GetMinDiagonal).
     * @param rA LHS matrix of the system.
     * @return Mean of the maximum and minimum diagonal values of rA.
     */
    static double GetAveragevalueDiagonal(const MatrixType& rA)
    {
        return 0.5 * (GetMaxDiagonal(rA) + GetMinDiagonal(rA));
    }

    /**
     * @brief This method returns the diagonal max value
     * @details Returns the diagonal max value of the LHS matrix rA, i.e. the largest absolute value of its diagonal
     * entries. It is evaluated in parallel by rows (a row without a stored diagonal entry does not contribute). The
     * matrix is first put in compressed mode if needed.
     * @param rA LHS matrix of the system.
     * @return Largest absolute value of the diagonal entries of rA.
     */
    static double GetMaxDiagonal(const MatrixType& rA)
    {
        EnsureCompressed(rA);
        const auto values = rA.value_data();
        const auto row_indices = rA.index1_data();
        const auto col_indices = rA.index2_data();

        return IndexPartition<std::size_t>(Size1(rA)).template for_each<MaxReduction<double>>([&](std::size_t Index) {
            const std::size_t col_begin = row_indices[Index];
            const std::size_t col_end = row_indices[Index + 1];
            for (std::size_t j = col_begin; j < col_end; ++j) {
                if (static_cast<std::size_t>(col_indices[j]) == Index) {
                    return std::abs(values[j]);
                }
            }
            return std::numeric_limits<double>::lowest();
        });
    }

    /**
     * @brief This method returns the diagonal min value
     * @details Returns the diagonal min value of the LHS matrix rA, i.e. the smallest absolute value of its
     * diagonal entries. It is evaluated in parallel by rows (a row without a stored diagonal entry does not
     * contribute). The matrix is first put in compressed mode if needed.
     * @param rA LHS matrix of the system.
     * @return Smallest absolute value of the diagonal entries of rA.
     */
    static double GetMinDiagonal(const MatrixType& rA)
    {
        EnsureCompressed(rA);
        const auto values = rA.value_data();
        const auto row_indices = rA.index1_data();
        const auto col_indices = rA.index2_data();

        return IndexPartition<std::size_t>(Size1(rA)).template for_each<MinReduction<double>>([&](std::size_t Index) {
            const std::size_t col_begin = row_indices[Index];
            const std::size_t col_end = row_indices[Index + 1];
            for (std::size_t j = col_begin; j < col_end; ++j) {
                if (static_cast<std::size_t>(col_indices[j]) == Index) {
                    return std::abs(values[j]);
                }
            }
            return std::numeric_limits<double>::max();
        });
    }

    ///@}
    ///@name Input and output
    ///@{

    /**
     * @brief Turn back information as a string
     * @details Returns the string "EigenSpace".
     * @return The string "EigenSpace".
     */
    virtual std::string Info() const
    {
        return "EigenSpace";
    }

    /**
     * @brief Print information about this object
     * @details Writes the string "EigenSpace" to the output stream rOStream.
     * @param rOStream Output stream where the information is written.
     */
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "EigenSpace";
    }

    /**
     * @brief Print object's data
     * @details Nothing is written to the output stream rOStream, as the space holds no data.
     * @param rOStream Output stream that would receive the data; nothing is written to it.
     */
    virtual void PrintData(std::ostream& rOStream) const
    {
    }

    //***********************************************************************

    /**
     * @brief Tells whether the space is distributed
     * @details It is a compile-time constant that always returns false, as this space is not distributed over MPI
     * ranks.
     * @return Always false, as the space is not distributed over MPI ranks.
     */
    inline static constexpr bool IsDistributed()
    {
        return false;
    }

    /**
     * @brief Returns a list of the fastest direct solvers
     * @details Returns a vector of strings with the names of the fastest direct solvers: pardiso_lu and
     * pardiso_ldlt (LinearSolversApplication, if compiled with Intel-support), sparse_lu (LinearSolversApplication)
     * and skyline_lu_factorization (in Core, always available, but slow). The order of the solvers in the list may
     * need to be updated and reordered depending on the size of the equation system.
     * @return Vector of strings with the names of the fastest direct solvers.
     */
    inline static std::vector<std::string> FastestDirectSolverList()
    {
        std::vector<std::string> faster_direct_solvers({
            "pardiso_lu",              // LinearSolversApplication (if compiled with Intel-support)
            "pardiso_ldlt",            // LinearSolversApplication (if compiled with Intel-support)
            "sparse_lu",               // LinearSolversApplication
            "skyline_lu_factorization" // In Core, always available, but slow
        });
        return faster_direct_solvers;
    }

    //***********************************************************************

    /**
     * @brief Returns an entry of a vector
     * @details Value of the entry I of the vector x.
     * @param x Vector from which the entry is read.
     * @param I Index of the entry.
     * @return Value of the entry I of x.
     */
    inline static TDataType GetValue(const VectorType& x, std::size_t I)
    {
        return x[I];
    }

    //***********************************************************************

    /**
     * @brief Gathers entries of a vector into an array
     * @details Copies the entries of the vector x at the positions given by IndexArray into the array pValues, so
     * that pValues[i] = x[IndexArray[i]]. The array pValues must have room for as many values as IndexArray has.
     * @param x Vector from which the entries are read.
     * @param IndexArray Positions of the entries of x to gather.
     * @param pValues Array that receives the gathered entries, pValues[i] = x[IndexArray[i]].
     */
    static void GatherValues(const VectorType& x, const std::vector<std::size_t>& IndexArray, TDataType* pValues)
    {
        KRATOS_TRY

        for (std::size_t i = 0; i < IndexArray.size(); i++)
            pValues[i] = x[IndexArray[i]];

        KRATOS_CATCH("")
    }

    /**
     * @brief Matrix-market output for the Eigen CSR matrix
     * @details The vendored Eigen has no unsupported/ module (no saveMarket), so the file is written directly from
     * the CSR arrays with the same format as matrix_market_interface.h. The matrix rM is written to the file
     * pFileName in coordinate format with one-based indices. If Symmetric is true, the symmetric banner is used and
     * only the entries on or below the diagonal are written. It returns true on success, and false if the file can
     * not be opened or the data can not be written.
     * @tparam TOtherMatrixType Type of the sparse matrix (an Eigen CSR matrix, e.g. EigenCompressedMatrix).
     * @param pFileName Name of the file to write.
     * @param rM Sparse matrix to write.
     * @param Symmetric If true, the symmetric banner is used and only the entries on or below the diagonal are written.
     * @return True on success, false if the file can not be opened or the data can not be written.
     */
    template<class TOtherMatrixType>
    static bool WriteMatrixMarketMatrix(const char* pFileName, const TOtherMatrixType& rM, const bool Symmetric)
    {
        FILE* f = fopen(pFileName, "w");
        if (f == NULL) {
            printf("WriteMatrixMarketMatrix(): unable to open %s.\n", pFileName);
            return false;
        }

        MM_typecode mm_code;
        mm_initialize_typecode(&mm_code);
        mm_set_matrix(&mm_code);
        mm_set_coordinate(&mm_code);
        SetMatrixMarketValueTypeCode(mm_code, typename TOtherMatrixType::value_type());
        if (Symmetric)
            mm_set_symmetric(&mm_code);
        else
            mm_set_general(&mm_code);
        mm_write_banner(f, mm_code);

        EnsureCompressed(rM);
        const auto* row_ptr = rM.outerIndexPtr();
        const auto* col_idx = rM.innerIndexPtr();
        const auto* values = rM.valuePtr();
        const std::size_t size1 = rM.size1();

        // Find out the actual number of non-zeros in case of a symmetric matrix
        int nnz = 0;
        if (Symmetric) {
            for (std::size_t i = 0; i < size1; ++i)
                for (auto k = row_ptr[i]; k < row_ptr[i + 1]; ++k)
                    if (i >= static_cast<std::size_t>(col_idx[k]))
                        nnz++;
        } else {
            nnz = static_cast<int>(rM.nnz());
        }

        mm_write_mtx_crd_size(f, rM.size1(), rM.size2(), nnz);

        for (std::size_t i = 0; i < size1; ++i) {
            for (auto k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
                const int I = static_cast<int>(i), J = static_cast<int>(col_idx[k]);
                if (Symmetric && I < J) continue;
                if (WriteMatrixMarketMatrixEntry(f, I + 1, J + 1, values[k]) < 0) {
                    printf("WriteMatrixMarketMatrix(): unable to write data.\n");
                    fclose(f);
                    return false;
                }
            }
        }

        fclose(f);
        return true;
    }

    /**
     * @brief Matrix-market output for a vector
     * @details Writes the vector rV to the file pFileName in the Matrix Market array format, as a column with
     * rV.size() rows. The value type of the banner (real or complex) follows the scalar type of the entries. It
     * returns true on success, and false if the file can not be opened or the data can not be written.
     * @tparam TOtherVectorType Type of the vector (it needs size and the call operator with an index).
     * @param pFileName Name of the file to write.
     * @param rV Vector to write.
     * @return True on success, false if the file can not be opened or the data can not be written.
     */
    template<class TOtherVectorType>
    static bool WriteMatrixMarketVector(const char* pFileName, const TOtherVectorType& rV)
    {
        FILE* f = fopen(pFileName, "w");
        if (f == NULL) {
            printf("WriteMatrixMarketVector(): unable to open %s.\n", pFileName);
            return false;
        }

        MM_typecode mm_code;
        mm_initialize_typecode(&mm_code);
        mm_set_matrix(&mm_code);
        mm_set_array(&mm_code);
        SetMatrixMarketValueTypeCode(mm_code, rV(0));
        mm_write_banner(f, mm_code);
        mm_write_mtx_array_size(f, rV.size(), 1);

        for (std::size_t i = 0; i < static_cast<std::size_t>(rV.size()); ++i) {
            if (WriteMatrixMarketVectorEntry(f, rV(i)) < 0) {
                printf("WriteMatrixMarketVector(): unable to write data.\n");
                fclose(f);
                return false;
            }
        }

        fclose(f);
        return true;
    }

    /**
     * @brief Creates a DoF updater
     * @details Returns a unique pointer to a new DoF updater for this space, created through DofUpdater::Create.
     * @return Unique pointer to the new DoF updater of this space.
     */
    static DofUpdaterPointerType CreateDofUpdater()
    {
        DofUpdaterType tmp;
        return tmp.Create();
    }

    ///@}
private:
    ///@name Private Operations
    ///@{

    /**
     * @brief Splits [0, Size) into one contiguous chunk per thread and calls f(Begin, End) in parallel
     * @details Keeping the hot loops in flat per-chunk bodies (instead of per-index lambdas inside the
     * OpenMP-outlined region) lets the compiler optimize them fully, matching the uBLAS performance. With a single
     * thread f is called once for the whole range, skipping the OpenMP region and the IndexPartition machinery,
     * whose fixed per-call cost is measurable on small and medium sizes.
     * @tparam TFunction Type of the callable invoked on each chunk as f(Begin, End).
     * @param Size Number of entries of the range [0, Size) to split.
     * @param f Callable invoked on each chunk with the bounds Begin and End (excluded) of the chunk.
     */
    template<class TFunction>
    static void ParallelChunks(const SizeType Size, TFunction&& f)
    {
        const int n_chunks = ParallelUtilities::GetNumThreads();
        // Serial short-circuit: skip the OpenMP region and the IndexPartition
        // machinery entirely (their fixed per-call cost is measurable on the
        // small/medium sizes when running single-threaded).
        if (n_chunks == 1) {
            f(SizeType(0), Size);
            return;
        }
        IndexPartition<int>(n_chunks, n_chunks).for_each([&](int Chunk) {
            f((Size * Chunk) / n_chunks, (Size * (Chunk + 1)) / n_chunks);
        });
    }

    /**
     * @brief Range portion of a dot product
     * @details It computes the non-conjugated dot product of the entries from Begin to End (excluded) of the arrays
     * pX and pY. The work is delegated to Eigen's redux over unaligned maps: it keeps several independent packet
     * accumulators, which beats both a strict scalar dependency chain and a manual 4-way unroll (the latter is
     * folded back into a single packet chain by the compiler). For complex scalars a transposed product is used, as
     * Eigen's dot() conjugates its first argument.
     * @param pX Array with the entries of the first vector.
     * @param pY Array with the entries of the second vector.
     * @param Begin Index of the first entry considered.
     * @param End Index of the entry after the last one considered (End is excluded).
     * @return Non-conjugated dot product of the entries from Begin to End (excluded) of pX and pY.
     */
    static TDataType PartialDot(
        const TDataType* pX,
        const TDataType* pY,
        const SizeType Begin,
        const SizeType End
        )
    {
        using MapType = Eigen::Map<const Eigen::Matrix<TDataType, Eigen::Dynamic, 1>, Eigen::Unaligned>;
        const auto n = static_cast<Eigen::Index>(End - Begin);
        if constexpr (Eigen::NumTraits<TDataType>::IsComplex) {
            // non-conjugated, as ublas inner_prod (Eigen's dot() conjugates its first argument)
            return (MapType(pX + Begin, n).transpose() * MapType(pY + Begin, n)).value();
        } else {
            return MapType(pX + Begin, n).dot(MapType(pY + Begin, n));
        }
    }

    /**
     * @brief Row-range portion of the sparse matrix-vector product
     * @details It is a plain function over raw pointers so the loop fully optimizes outside the OpenMP-outlined
     * caller. For each row i from Begin to End (excluded), pY[i] is set to the sum of pValues[k] times
     * pX[pColumnIndices[k]] over the entries k of the row given by the CSR row pointers pRowIndices.
     * @tparam TIndexType Index type of the CSR row pointers and column indices.
     * @param pRowIndices CSR row pointers of the matrix.
     * @param pColumnIndices CSR column indices of the matrix.
     * @param pValues CSR values of the matrix.
     * @param pX Array with the entries of the vector multiplied by the matrix.
     * @param pY Array that receives the entries of the product for the rows from Begin to End (excluded).
     * @param Begin Index of the first row computed.
     * @param End Index of the row after the last one computed (End is excluded).
     */
    template<class TIndexType>
    static void PartialSpMV(
        const TIndexType* pRowIndices,
        const TIndexType* pColumnIndices,
        const TDataType* pValues,
        const TDataType* pX,
        TDataType* pY,
        const SizeType Begin,
        const SizeType End
        )
    {
        for (SizeType i = Begin; i < End; ++i) {
            TDataType accumulator = TDataType();
            const auto row_end = pRowIndices[i + 1];
            for (auto k = pRowIndices[i]; k < row_end; ++k) {
                accumulator += pValues[k] * pX[pColumnIndices[k]];
            }
            pY[i] = accumulator;
        }
    }

    /**
     * @brief Makes sure that a sparse matrix is in compressed mode
     * @details The CSR-array based operations require the matrix to be in compressed mode. A plain Eigen sparse matrix
     * can be left uncompressed by insertion, and it is then normalized here (logically const: the represented matrix
     * does not change). The operator() of EigenCompressedMatrix keeps its storage packed, so it does not normally need
     * this call.
     * @tparam TOtherMatrixType Type of the sparse matrix (it needs isCompressed and makeCompressed).
     * @param rA Sparse matrix to put in compressed mode if needed.
     */
    template<class TOtherMatrixType>
    static void EnsureCompressed(const TOtherMatrixType& rA)
    {
        if (!rA.isCompressed()) {
            const_cast<TOtherMatrixType&>(rA).makeCompressed();
        }
    }

    ///@}
    ///@name Un accessible methods
    ///@{

    /**
     * @brief Assignment operator
     * @details It is private and not implemented, so the space can not be assigned.
     * @param rOther Space to assign from.
     * @return Reference to this space (the operator is not implemented).
     */
    EigenSpace& operator=(EigenSpace const& rOther);

    /**
     * @brief Copy constructor
     * @details It is private and not implemented, so the space can not be copied.
     * @param rOther Space to copy.
     */
    EigenSpace(EigenSpace const& rOther);

    ///@}
}; // Class EigenSpace

///@}

} // namespace Kratos.
