//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#pragma once

// System includes

// External includes

// Project includes
#include "containers/csr_matrix.h"
#include "includes/model_part.h"

namespace Kratos::Future
{

///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

namespace
{
    // Lazy wrapper to prevent instantiation unless selected
    template<typename... Ts>
    struct matrix_type_or_void
    {
        using type = void;
        using p_type = void*;
    };

    template<typename T0, typename... Ts>
    struct matrix_type_or_void<T0, Ts...>
    {
        using type = T0;
        using p_type = typename T0::Pointer;
    };
}

/**
 * @brief
 * @class LinearOperator
 * @ingroup KratosCore
 * @brief Auxiliary container to store a linear operator
 * @details
 * This lightweight container mimics the behavior of a linear operator,
 * providing a general abstraction for the application of a matrix-like object
 * onto a vector without explicitly storing the matrix entries.
 *
 * It is conceptually equivalent to the Python `scipy.sparse.linalg.LinearOperator`.
 * The operator defines two callable functions:
 *  - `Apply()`: performs the product \( y = A(x) \)
 *  - `ApplyTranspose()`: performs the product \( y = A^T(x) \)
 *
 * This structure can be used to:
 *  - Represent explicit CSR matrices (via wrapping)
 *  - Derive matrix-free operators (e.g. preconditioners, Schur complements)
 * @tparam TSystemVectorType The system vector type used in the Kratos linear algebra backend
 * @tparam Ts Optional parameter pack representing the sparse matrix type (e.g., CsrMatrix<double>).
 */
template <class TSystemVectorType, typename... Ts>
class LinearOperator
{

static_assert(sizeof...(Ts) < 2, "Ts must be either a single type representing the sparse matrix type or empty to indicate a matrix-free operator.");

public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of LinearOperator
    KRATOS_CLASS_POINTER_DEFINITION(LinearOperator);

    /// Sparse matrix type (if any) or void
    using TSparseMatrixType = typename matrix_type_or_void<Ts...>::type;

    /// Sparse matrix pointer type (if any) or void*
    using TSparseMatrixTypePointer = typename matrix_type_or_void<Ts...>::p_type;

    /// Boolean indicating whether a CSR matrix type is provided
    static constexpr bool kIsMatrixFree = std::is_void_v<TSparseMatrixType>;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @details Creates an empty LinearOperator with uninitialized function objects.
     */
    LinearOperator() = default;

    /**
     * @brief Constructor from row/column sizes.
     * @param NumRows Number of rows of the operator
     * @param NumCols Number of columns of the operator
     */
    LinearOperator(
        std::size_t NumRows,
        std::size_t NumCols) requires(kIsMatrixFree)
        : mNumRows(NumRows)
        , mNumCols(NumCols)
    {
    }

    /**
     * @brief Constructor from a pair of (rows, cols).
     * @param Size Pair containing {NumRows, NumCols}
     */
    LinearOperator(std::pair<std::size_t, std::size_t> Size) requires(kIsMatrixFree)
        : mNumRows(std::get<0>(Size))
        , mNumCols(std::get<1>(Size))
    {
    }

    /**
     * @brief Constructor from a CSR matrix.
     * Constructs a LinearOperator that wraps the provided CSR matrix.
     * Note that this constructor is only enabled when a CSR matrix type is provided.
     * @param rA Reference to the CSR matrix
     */
    LinearOperator(Ts& ...rA) requires(!kIsMatrixFree)
    {
        auto& r_A = std::get<0>(std::forward_as_tuple(rA...)); // Extract the CSR matrix from the parameter pack
        mpCsrMatrix = std::make_shared<TSparseMatrixType>(r_A);
        mNumRows = mpCsrMatrix->size1();
        mNumCols = mpCsrMatrix->size2();
    }

    /// Deleted copy constructor (non-copyable)
    LinearOperator(const LinearOperator& rOther) = delete;

    /// Defaulted move constructor
    LinearOperator(LinearOperator&& rOther) = default;

    /// Default destructor (non-virtual)
    ~LinearOperator() = default;

    ///@}
    ///@name Operators
    ///@{

    /// Deleted copy assignment operator (non-copyable)
    LinearOperator& operator=(const LinearOperator& rOther) = delete;

    /// Defaulted move assignment operator
    LinearOperator& operator=(LinearOperator&& rOther) = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Apply the linear operator.
     * @details Computes \( y = A(x) \).
     * @param rX Input vector
     * @param rY Output vector (result)
     */
    void Apply(
        const TSystemVectorType& rX,
        TSystemVectorType& rY)
    {
        if constexpr (!kIsMatrixFree) {
            // Apply the CSR matrix
            mpCsrMatrix->SpMV(rX, rY);
        } else {
            // Matrix-free application (to be implemented by the user)
            KRATOS_ERROR << "Matrix-free Apply() must be implemented in derived classes." << std::endl;
        }
    }

    /**
     * @brief Apply the transpose of the linear operator.
     * @details Computes \( y = A^T(x) \).
     * @param rX Input vector
     * @param rY Output vector (result)
     */
    void ApplyTranspose(
        const TSystemVectorType& rX,
        TSystemVectorType& rY)
    {
        if constexpr (!kIsMatrixFree) {
            // Apply the transpose of the CSR matrix
            mpCsrMatrix->TransposeSpMV(rX, rY);
        } else {
            // Matrix-free transpose application (to be implemented by the user)
            KRATOS_ERROR << "Matrix-free ApplyTranspose() must be implemented in derived classes." << std::endl;
        }
    }

    /**
     * @brief Clear the operator data.
     * @details Resets the sizes and function objects to null.
     */
    void Clear()
    {
        mNumRows = 0;
        mNumCols = 0;
        mpCsrMatrix = nullptr;
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Set the number of rows.
     * @param NumRows Number of rows
     */
    void SetRows(std::size_t NumRows)
    {
        mNumRows = NumRows;
    }

    /**
     * @brief Set the number of columns.
     * @param NumCols Number of columns
     */
    void SetCols(std::size_t NumCols)
    {
        mNumCols = NumCols;
    }

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief Get the number of rows.
     * @return Number of rows of the operator
     */
    std::size_t Rows() const
    {
        return mNumRows;
    }

    /**
     * @brief Get the number of columns.
     * @return Number of columns of the operator
     */
    std::size_t Cols() const
    {
        return mNumCols;
    }

    /**
     * @brief Check if the operator is matrix-free.
     * @return True if the operator is matrix-free, false if it wraps a CSR matrix
     */
    bool IsMatrixFree() const
    {
        return kIsMatrixFree;
    }

    ///@}

private:

    /// Number of rows of the operator
    std::size_t mNumRows = 0;

    /// Number of columns of the operator
    std::size_t mNumCols = 0;

    /// Pointer to the CSR matrix (if applicable)
    TSparseMatrixTypePointer mpCsrMatrix = nullptr;

}; // class LinearOperator

///@}
///@name Input and output
///@{

///@}
///@} addtogroup block

} // namespace Kratos::Future
