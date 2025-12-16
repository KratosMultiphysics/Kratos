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
#include "containers/distributed_csr_matrix.h"
#include "future/containers/linear_operator.h"
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

// /**
//  * @brief
//  * @class SparseMatrixLinearOperator
//  * @ingroup KratosCore
//  * @brief Auxiliary container to store a linear operator
//  * @details
//  * This lightweight container mimics the behavior of a linear operator,
//  * providing a general abstraction for the application of a matrix-like object
//  * onto a vector without explicitly storing the matrix entries.
//  *
//  * It is conceptually equivalent to the Python `scipy.sparse.linalg.SparseMatrixLinearOperator`.
//  * The operator defines two callable functions:
//  *  - `Apply()`: performs the product \( y = A(x) \)
//  *  - `ApplyTranspose()`: performs the product \( y = A^T(x) \)
//  *
//  * This structure can be used to:
//  *  - Represent explicit CSR matrices (via wrapping)
//  *  - Derive matrix-free operators (e.g. preconditioners, Schur complements)
//  * @tparam TVectorType The system vector type used in the Kratos linear algebra backend
//  * @tparam Ts Optional parameter pack representing the sparse matrix type (e.g., CsrMatrix<double>).
//  */
template <class TVectorType>
class SparseMatrixLinearOperator : public LinearOperator<TVectorType>
{

public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SparseMatrixLinearOperator
    KRATOS_CLASS_POINTER_DEFINITION(SparseMatrixLinearOperator);

    /// Data type stored in the system vector
    using DataType = typename TVectorType::DataType;

    /// Index type used in the system vector
    using IndexType = typename TVectorType::IndexType;

    /// Serial CSR matrix type
    using CsrMatrixType = CsrMatrix<DataType, IndexType>;

    /// Distributed CSR matrix type
    using DistributedCsrMatrixType = DistributedCsrMatrix<DataType, IndexType>;

    /// Matrix variant type definition from base class
    using MatrixVariantType = typename LinearOperator<TVectorType>::MatrixVariantType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @details Creates an empty SparseMatrixLinearOperator with uninitialized function objects.
     */
    SparseMatrixLinearOperator() = default;

    /**
     * @brief Constructor from parameters.
     * @param ThisParameters Parameters containing the linear operator settings
     */
    SparseMatrixLinearOperator(Parameters ThisParameters)
        : LinearOperator<TVectorType>(ThisParameters)
    {
    }

    /**
     * @brief Constructor from a CSR matrix.
     * Constructs a SparseMatrixLinearOperator that wraps the provided CSR matrix.
     * Note that this constructor is only enabled when a CSR matrix type is provided.
     * @param rA Reference to the CSR matrix
     */
    template<class TMatrixType>
    SparseMatrixLinearOperator(const TMatrixType& rA)
        : LinearOperator<TVectorType>()
    {
        mCsrMatrix = rA;
        this->SetNumRows(rA.size1());
        this->SetNumCols(rA.size2());
    }

    /// Deleted copy constructor (non-copyable)
    SparseMatrixLinearOperator(const SparseMatrixLinearOperator& rOther) = delete;

    /// Defaulted move constructor
    SparseMatrixLinearOperator(SparseMatrixLinearOperator&& rOther) = default;

    /// Default destructor (non-virtual)
    ~SparseMatrixLinearOperator() = default;

    ///@}
    ///@name Operators
    ///@{

    /// Deleted copy assignment operator (non-copyable)
    SparseMatrixLinearOperator& operator=(const SparseMatrixLinearOperator& rOther) = delete;

    /// Defaulted move assignment operator
    SparseMatrixLinearOperator& operator=(SparseMatrixLinearOperator&& rOther) = default;

    ///@}
    ///@name Operations
    ///@{

    void SpMV(
        const TVectorType& rX,
        TVectorType& rY) override
    {
        if (std::holds_alternative<>(mCsrMatrix)) {
            std::get<CsrMatrixType>(mCsrMatrix).SpMV(rX, rY);
        } else if (std::holds_alternative<DistributedCsrMatrixType>(mCsrMatrix)) {
            std::get<DistributedCsrMatrixType>(mCsrMatrix).SpMV(rX, rY);
        } else {
            KRATOS_ERROR << "Unsupported matrix type in SparseMatrixLinearOperator Apply()." << std::endl;
        }
    }

    void TransposeSpMV(
        const TVectorType& rX,
        TVectorType& rY) override
    {
        if (std::holds_alternative<CsrMatrixType>(mCsrMatrix)) {
            std::get<CsrMatrixType>(mCsrMatrix).TransposeSpMV(rX, rY);
        } else if (std::holds_alternative<DistributedCsrMatrixType>(mCsrMatrix)) {
            std::get<DistributedCsrMatrixType>(mCsrMatrix).TransposeSpMV(rX, rY);
        } else {
            KRATOS_ERROR << "Unsupported matrix type in SparseMatrixLinearOperator Apply()." << std::endl;
        }
    }

    void Clear() override
    {
        LinearOperator<TVectorType>::Clear();
        mCsrMatrix = std::monostate{};
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    bool IsMatrixFree() const override
    {
        return false;
    }

    ///@}

protected:

    ///@name Protected access
    ///@{

    MatrixVariantType& GetMatrixImpl() const override
    {
        return mCsrMatrix;
    }

    ///@}

private:

    /// Pointer to the CSR matrix (if applicable)
    MatrixVariantType mCsrMatrix;

}; // class SparseMatrixLinearOperator

///@}
///@name Input and output
///@{

///@}
///@} addtogroup block

} // namespace Kratos::Future
