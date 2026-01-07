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
#include <any>

// External includes

// Project includes
#include "containers/csr_matrix.h"
#include "containers/distributed_csr_matrix.h"
#include "includes/model_part.h"

namespace Kratos::Future
{

///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

// namespace
// {
//     // Lazy wrapper to prevent instantiation unless selected
//     template<typename... Ts>
//     struct matrix_type_or_void
//     {
//         using type = void;
//         using p_type = void*;
//     };

//     template<typename T0, typename... Ts>
//     struct matrix_type_or_void<T0, Ts...>
//     {
//         using type = T0;
//         using p_type = typename T0::Pointer;
//     };
// }

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
template <class TVectorType>
class LinearOperator
{

public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of LinearOperator
    KRATOS_CLASS_POINTER_DEFINITION(LinearOperator);

    /// Data type stored in the system vector
    using DataType = typename TVectorType::DataType;

    /// Index type used in the system vector
    using IndexType = typename TVectorType::IndexType;

    /// Any type used to represent system matrix
    using MatrixType = std::any;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @details Creates an empty LinearOperator with uninitialized function objects.
     */
    LinearOperator() = default;

    /**
     * @brief Constructor from parameters.
     * @param ThisParameters Parameters containing the linear operator settings
     */
    LinearOperator(Parameters ThisParameters)
    {
        mNumRows = ThisParameters["num_rows"].GetInt();
        mNumCols = ThisParameters["num_cols"].GetInt();
    }

    /// Deleted copy constructor (non-copyable)
    LinearOperator(const LinearOperator& rOther) = delete;

    /// Defaulted move constructor
    LinearOperator(LinearOperator&& rOther) = default;

    /// Default destructor (non-virtual)
    virtual ~LinearOperator() = default;

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

    virtual void SpMV(
        const TVectorType& rX,
        TVectorType& rY) const = 0;

    virtual void TransposeSpMV(
        const TVectorType& rX,
        TVectorType& rY) const = 0;

    /**
     * @brief Clear the operator data.
     * @details Resets the sizes and function objects to null.
     */
    virtual void Clear()
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
    void SetNumRows(std::size_t NumRows)
    {
        mNumRows = NumRows;
    }

    /**
     * @brief Set the number of columns.
     * @param NumCols Number of columns
     */
    void SetNumCols(std::size_t NumCols)
    {
        mNumCols = NumCols;
    }

    template<class TMatrixType>
    TMatrixType& GetMatrix()
    {
        auto& r_matrix = this->GetMatrixImpl();
        // KRATOS_ERROR_IF_NOT(std::holds_alternative<TMatrixType>(r_matrix)) << "Requested matrix type does not match stored matrix.";
        // return std::get<TMatrixType>(r_matrix);
        return std::any_cast<TMatrixType&>(r_matrix);
    }

    template<class TMatrixType>
    const TMatrixType& GetMatrix() const
    {
        const auto& r_matrix = this->GetMatrixImpl();
        // KRATOS_ERROR_IF_NOT(std::holds_alternative<TMatrixType>(r_matrix)) << "Requested matrix type does not match stored matrix.";
        // return std::get<TMatrixType>(r_matrix);
        return std::any_cast<const TMatrixType&>(r_matrix);
    }

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief Get the number of rows.
     * @return Number of rows of the operator
     */
    std::size_t NumRows() const //TODO: Rename to size1() as in our CSR arrays?
    {
        return mNumRows;
    }

    /**
     * @brief Get the number of columns.
     * @return Number of columns of the operator
     */
    std::size_t NumCols() const //TODO: Rename to size2() as in our CSR arrays?
    {
        return mNumCols;
    }

    /**
     * @brief Check if the operator is matrix-free.
     * @return True if the operator is matrix-free, false if it wraps a CSR matrix
     */
    virtual bool IsMatrixFree() const
    {
        return true;
    }

    ///@}

protected:

    ///@name Protected access
    ///@{

    /**
     * @brief Get the underlying matrix variant.
     * @return Reference to the matrix variant
     */
    [[noreturn]] virtual MatrixType& GetMatrixImpl()
    {
        KRATOS_ERROR << "GetMatrixImpl() not implemented in base LinearOperator class." << std::endl;
    }
    /**
     * @brief Get the underlying matrix variant.
     * @return Reference to the matrix variant
     */
    [[noreturn]] virtual const MatrixType& GetMatrixImpl() const
    {
        KRATOS_ERROR << "GetMatrixImpl() not implemented in base LinearOperator class." << std::endl;
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
