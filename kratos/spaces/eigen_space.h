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
#include "includes/kratos_eigen_interface.h"
#include "includes/matrix_market_interface.h"
#include "spaces/ublas_space.h" // For the SCALING_DIAGONAL enum and the dense Kratos::Matrix used in Get/SetColumn
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

template <class TDataType, class TMatrixType, class TVectorType>
class EigenSpace;

template <class TDataType, class TIndexType = KratosEigenIndexType>
using TEigenSparseSpace =
    EigenSpace<TDataType, EigenCompressedMatrix<TDataType, TIndexType>, EigenVector<TDataType>>;

template <class TDataType>
using TEigenDenseSpace =
    EigenSpace<TDataType, EigenMatrix<TDataType>, EigenVector<TDataType>>;

///@}
///@name Kratos Classes
///@{

/**
 * @class EigenSpace
 * @ingroup KratosCore
 * @brief A class template for handling data types, matrices, and vectors in an Eigen-based space.
 * @details Mirrors the static API of UblasSpace over the Eigen backend types
 * (EigenCompressedMatrix, EigenMatrix, EigenVector), so it can be used as
 * TSparseSpace/TDenseSpace template argument by the strategies, builder-and-solvers
 * and linear solvers.
 * @tparam TDataType The data type used in the Eigen space.
 * @tparam TMatrixType The matrix type used in the Eigen space.
 * @tparam TVectorType The vector type used in the Eigen space.
 * @author Vicente Mataix Ferrandiz
 */
template<class TDataType, class TMatrixType, class TVectorType>
class EigenSpace
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of EigenSpace
    KRATOS_CLASS_POINTER_DEFINITION(EigenSpace);

    /// The data type considered
    using DataType = TDataType;

    /// The matrix type considered
    using MatrixType = TMatrixType;

    /// The vector type considered
    using VectorType = TVectorType;

    /// The index type considered
    using IndexType = std::size_t;

    /// The size type considered
    using SizeType = std::size_t;

    /// The pointer to the matrix type
    using MatrixPointerType = typename Kratos::shared_ptr<TMatrixType>;

    /// The pointer to the vector type
    using VectorPointerType = typename Kratos::shared_ptr<TVectorType>;

    /// The DoF updater type
    using DofUpdaterType = DofUpdater<EigenSpace<TDataType, TMatrixType, TVectorType>>;

    /// The pointer to the DoF updater type
    using DofUpdaterPointerType = typename DofUpdaterType::UniquePointer;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    EigenSpace()
    {
    }

    /// Destructor.
    virtual ~EigenSpace()
    {
    }

    ///@}
    ///@name Operations
    ///@{

    static MatrixPointerType CreateEmptyMatrixPointer()
    {
        return MatrixPointerType(new TMatrixType(0, 0));
    }

    static VectorPointerType CreateEmptyVectorPointer()
    {
        return VectorPointerType(new TVectorType(0));
    }

    /**
     * @brief This method returns true if the pointer is null
     * @param pPointer The pointer considered
     * @return True if the pointer is null
     */
    template<class TPointerType>
    inline static bool IsNull(const TPointerType& pPointer)
    {
        return !static_cast<bool>(pPointer);
    }

    /// return size of vector rV
    static IndexType Size(VectorType const& rV)
    {
        return static_cast<IndexType>(rV.size());
    }

    /// return number of rows of rM
    static IndexType Size1(MatrixType const& rM)
    {
        return rM.size1();
    }

    /// return number of columns of rM
    static IndexType Size2(MatrixType const& rM)
    {
        return rM.size2();
    }

    /// rXi = rMij (works on the dense Kratos::Matrix used for multi-column solves, as in UblasSpace)
    template<typename TColumnType>
    static void GetColumn(unsigned int j, Matrix& rM, TColumnType& rX)
    {
        if (rX.size() != rM.size1())
            rX.resize(rM.size1(), false);

        for (IndexType i = 0; i < rM.size1(); i++) {
            rX[i] = rM(i, j);
        }
    }

    template<typename TColumnType>
    static void SetColumn(unsigned int j, Matrix& rM, TColumnType& rX)
    {
        for (IndexType i = 0; i < rM.size1(); i++) {
            rM(i, j) = rX[i];
        }
    }

    /// rY = rX
    static void Copy(MatrixType const& rX, MatrixType& rY)
    {
        rY = rX;
    }

    /// rY = rX
    static void Copy(VectorType const& rX, VectorType& rY)
    {
        const SizeType size = Size(rX);
        if (Size(rY) != size)
            rY.resize(size, false);
        const auto* x_data = rX.data();
        auto* y_data = rY.data();
        ParallelChunks(size, [=](const SizeType Begin, const SizeType End) {
            for (SizeType i = Begin; i < End; ++i)
                y_data[i] = x_data[i];
        });
    }

    /// rX * rY (non-conjugated, as in UblasSpace)
    static TDataType Dot(VectorType const& rX, VectorType const& rY)
    {
        const SizeType size = Size(rX);
        const auto* x_data = rX.data();
        const auto* y_data = rY.data();
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

    /// ||rX||2
    static TDataType TwoNorm(VectorType const& rX)
    {
        return std::sqrt(Dot(rX, rX));
    }

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

    /// TwoNorm of a (dense) Eigen expression, e.g. TwoNorm(y - x). The exact
    /// template match avoids the ambiguity between the vector and dense-matrix
    /// overloads that a converting expression argument would otherwise cause.
    template<class TDerived>
    static TDataType TwoNorm(const Eigen::MatrixBase<TDerived>& rX)
    {
        return rX.norm();
    }

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
     * @param rA The matrix to compute the Jacobi norm
     * @return The Jacobi norm
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

    static void Mult(const EigenMatrix<TDataType>& rA, const VectorType& rX, VectorType& rY)
    {
        if (Size(rY) != rA.size1())
            rY.resize(rA.size1(), false);
        rY.noalias() = rA * rX;
    }

    /// Sparse SpMV. Deterministic manual CSR loop (as in KratosSpace) rather than
    /// Eigen's internally parallelized product, to avoid nested parallelism under
    /// the Kratos shared-memory parallelization and keep UblasSpace behavior.
    /// The per-chunk work is delegated to a plain function taking the raw
    /// pointers as arguments (as UblasSpace does): keeping the hot loop out of
    /// the OpenMP-outlined lambda lets the compiler optimize it fully, which
    /// benchmarks at parity with the uBLAS implementation.
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
        const auto* x_data = rX.data();
        auto* y_data = rY.data();
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
    template<class TOtherMatrixType>
    static void TransposeMult(const TOtherMatrixType& rA, const VectorType& rX, VectorType& rY)
    {
        if (Size(rY) != rA.size2())
            rY.resize(rA.size2(), false);
        rY.noalias() = rA.transpose() * rX;
    }

    /// Sparse transpose product: a deterministic serial CSR scatter
    /// (y[col] += v * x[row]), exactly the access pattern of the uBLAS
    /// axpy_prod(x, A, y) implementation (which is serial too). Eigen's
    /// expression path for this case does an extra full setZero() pass over rY
    /// and re-tests isCompressed() per row inside its InnerIterator, which
    /// benchmarks measurably slower.
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
        const auto* x_data = rX.data();
        auto* y_data = rY.data();

        SetToZero(rY);
        for (SizeType i = 0; i < n_rows; ++i) {
            const TDataType x_i = x_data[i];
            const auto row_end = row_ptr[i + 1];
            for (auto k = row_ptr[i]; k < row_end; ++k) {
                y_data[col_idx[k]] += values[k] * x_i;
            }
        }
    }

    static inline SizeType GraphDegree(IndexType i, TMatrixType& rA)
    {
        EnsureCompressed(rA);
        const auto* row_ptr = rA.outerIndexPtr();
        return static_cast<SizeType>(row_ptr[i + 1] - row_ptr[i]);
    }

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
    static void InplaceMult(VectorType& rX, const double A)
    {
        if (A == 1.00)
            return;
        auto* x_data = rX.data();
        ParallelChunks(Size(rX), [=](const SizeType Begin, const SizeType End) {
            for (SizeType i = Begin; i < End; ++i)
                x_data[i] *= A;
        });
    }

    //********************************************************************
    // ATTENTION it is assumed no aliasing between rX and rY
    // X = A*y;
    static void Assign(VectorType& rX, const double A, const VectorType& rY)
    {
        const SizeType size = Size(rY);
        if (Size(rX) != size)
            rX.resize(size, false);

        const auto* y_data = rY.data();
        auto* x_data = rX.data();
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
    static void UnaliasedAdd(VectorType& rX, const double A, const VectorType& rY)
    {
        const SizeType size = Size(rY);
        if (Size(rX) != size)
            rX.resize(size, false);

        const auto* y_data = rY.data();
        auto* x_data = rX.data();
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

    static void ScaleAndAdd(const double A, const VectorType& rX, const double B, const VectorType& rY, VectorType& rZ) // rZ = (A * rX) + (B * rY)
    {
        const SizeType size = Size(rX);
        if (Size(rZ) != size)
            rZ.resize(size, false);

        const auto* x_data = rX.data();
        const auto* y_data = rY.data();
        auto* z_data = rZ.data();
        ParallelChunks(size, [=](const SizeType Begin, const SizeType End) {
            for (SizeType i = Begin; i < End; ++i)
                z_data[i] = A * x_data[i] + B * y_data[i];
        });
    }

    static void ScaleAndAdd(const double A, const VectorType& rX, const double B, VectorType& rY) // rY = (A * rX) + (B * rY)
    {
        const SizeType size = Size(rX);
        const auto* x_data = rX.data();
        auto* y_data = rY.data();
        ParallelChunks(size, [=](const SizeType Begin, const SizeType End) {
            for (SizeType i = Begin; i < End; ++i)
                y_data[i] = A * x_data[i] + B * y_data[i];
        });
    }

    /// rA[i] * rX (non-conjugated, as in UblasSpace)
    static double RowDot(unsigned int i, MatrixType& rA, VectorType& rX)
    {
        TDataType accumulator = TDataType();
        for (IndexType j = 0; j < rA.size2(); ++j) {
            accumulator += rA(i, j) * rX[j];
        }
        return accumulator;
    }

    static void SetValue(VectorType& rX, IndexType i, TDataType value)
    {
        rX[i] = value;
    }

    /// rX = A
    static void Set(VectorType& rX, TDataType A)
    {
        rX.setConstant(A);
    }

    static void Resize(MatrixType& rA, SizeType m, SizeType n)
    {
        rA.resize(m, n, false);
    }

    static void Resize(MatrixPointerType& pA, SizeType m, SizeType n)
    {
        pA->resize(m, n, false);
    }

    static void Resize(VectorType& rX, SizeType n)
    {
        rX.resize(n, false);
    }

    static void Resize(VectorPointerType& pX, SizeType n)
    {
        pX->resize(n, false);
    }

    static void Clear(MatrixPointerType& pA)
    {
        pA->resize(0, 0, false);
    }

    static void Clear(VectorPointerType& pX)
    {
        pX->resize(0, false);
    }

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

    inline static void ResizeData(VectorType& rX, SizeType m)
    {
        rX.resize(m, false);
        rX.setZero();
    }

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

    inline static void SetToZero(VectorType& rX)
    {
        auto* x_data = rX.data();
        ParallelChunks(Size(rX), [=](const SizeType Begin, const SizeType End) {
            for (SizeType i = Begin; i < End; ++i)
                x_data[i] = TDataType();
        });
    }

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
     * @details This method returns the scale norm considering for scaling the diagonal
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
     * @param rA The LHS matrix
     * @return The diagonal norm
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
     * @param rA The LHS matrix
     * @return The diagonal average value
     */
    static double GetAveragevalueDiagonal(const MatrixType& rA)
    {
        return 0.5 * (GetMaxDiagonal(rA) + GetMinDiagonal(rA));
    }

    /**
     * @brief This method returns the diagonal max value
     * @param rA The LHS matrix
     * @return The diagonal max value
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
     * @param rA The LHS matrix
     * @return The diagonal min value
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

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "EigenSpace";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "EigenSpace";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }

    //***********************************************************************

    inline static constexpr bool IsDistributed()
    {
        return false;
    }

    /**
     * @brief Returns a list of the fastest direct solvers.
     * @details This function returns a vector of strings representing the names of the fastest direct solvers. The order of the solvers in the list may need to be updated and reordered depending on the size of the equation system.
     * @return A vector of strings containing the names of the fastest direct solvers.
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

    inline static TDataType GetValue(const VectorType& x, std::size_t I)
    {
        return x[I];
    }

    //***********************************************************************

    static void GatherValues(const VectorType& x, const std::vector<std::size_t>& IndexArray, TDataType* pValues)
    {
        KRATOS_TRY

        for (std::size_t i = 0; i < IndexArray.size(); i++)
            pValues[i] = x[IndexArray[i]];

        KRATOS_CATCH("")
    }

    /// Matrix-market output for the Eigen CSR matrix. The vendored Eigen has
    /// no unsupported/ module (no saveMarket), so the file is written directly
    /// from the CSR arrays with the same format as matrix_market_interface.h.
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

    static DofUpdaterPointerType CreateDofUpdater()
    {
        DofUpdaterType tmp;
        return tmp.Create();
    }

    ///@}
private:
    ///@name Private Operations
    ///@{

    /// Splits [0, Size) into one contiguous chunk per thread and calls
    /// f(Begin, End) in parallel. Keeping the hot loops in flat per-chunk
    /// bodies (instead of per-index lambdas inside the OpenMP-outlined region)
    /// lets the compiler optimize them fully, matching the uBLAS performance.
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

    /// Range portion of a dot product, delegated to Eigen's redux over
    /// unaligned maps: it keeps several independent packet accumulators, which
    /// beats both a strict scalar dependency chain and a manual 4-way unroll
    /// (the latter is folded back into a single packet chain by the compiler).
    static TDataType PartialDot(
        const TDataType* pX,
        const TDataType* pY,
        const SizeType Begin,
        const SizeType End
        )
    {
        using MapType = Eigen::Map<const Eigen::Matrix<TDataType, Eigen::Dynamic, 1>, Eigen::Unaligned>;
        const auto n = static_cast<Eigen::Index>(End - Begin);
        return MapType(pX + Begin, n).dot(MapType(pY + Begin, n));
    }

    /// Row-range portion of the sparse matrix-vector product; a plain function
    /// over raw pointers so the loop fully optimizes outside the OpenMP-outlined
    /// caller.
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

    /// The CSR-array based operations require the matrix to be in compressed
    /// mode. Element insertion through operator() (as done e.g. from python or
    /// by AssembleLHS) may leave it uncompressed, so it is normalized here
    /// (logically const: the represented matrix does not change).
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

    /// Assignment operator.
    EigenSpace& operator=(EigenSpace const& rOther);

    /// Copy constructor.
    EigenSpace(EigenSpace const& rOther);

    ///@}
}; // Class EigenSpace

///@}

} // namespace Kratos.
