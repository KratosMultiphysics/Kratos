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

/**
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
 * This structure can be used to represent:
 *  - Explicit matrices (via wrapping)
 *  - Matrix-free operators (e.g. preconditioners, Schur complements)
 *  - Transpose operations
 *
 * @tparam TSystemVectorType The system vector type used in the Kratos linear algebra backend
 */
template <class TSystemVectorType, typename... TSparseMatrixType>
class LinearOperator
{

static_assert(sizeof...(TSparseMatrixType) < 2, "TSparseMatrixType must be either a single type representing the sparse matrix type or empty to indicate a matrix-free operator.");

public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of LinearOperator
    KRATOS_CLASS_POINTER_DEFINITION(LinearOperator);

    static constexpr bool IsCsr = sizeof...(TSparseMatrixType) == 1;

    using TSparseMatrixTypePointer = std::conditional_t<IsCsr, typename std::tuple_element<0, std::tuple<TSparseMatrixType...>>::type::Pointer, void*>;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @details Creates an empty LinearOperator with uninitialized function objects.
     */
    LinearOperator() = default;

    /**
     * @brief Constructor from row/column sizes and function objects.
     * @param NumRows Number of rows of the operator
     * @param NumCols Number of columns of the operator
     */
    LinearOperator(
        std::size_t NumRows,
        std::size_t NumCols) requires(!IsCsr)
        : mNumRows(NumRows)
        , mNumCols(NumCols)
    {
    }

    /**
     * @brief Constructor from a pair of (rows, cols) and function objects.
     * @param Size Pair containing {NumRows, NumCols}
     */
    LinearOperator(std::pair<std::size_t, std::size_t> Size) requires(!IsCsr)
        : mNumRows(std::get<0>(Size))
        , mNumCols(std::get<1>(Size))
    {
    }

    LinearOperator(TSparseMatrixType... Args) requires(IsCsr)
    {
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
    }

    /**
     * @brief Clear the operator data.
     * @details Resets the sizes and function objects to null.
     */
    void Clear()
    {
        mNumRows = 0;
        mNumCols = 0;
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

    ///@}

private:

    /// Number of rows of the operator
    std::size_t mNumRows = 0;

    /// Number of columns of the operator
    std::size_t mNumCols = 0;

}; // class LinearOperator

///@}
///@name Input and output
///@{

///@}
///@} addtogroup block

} // namespace Kratos::Future
