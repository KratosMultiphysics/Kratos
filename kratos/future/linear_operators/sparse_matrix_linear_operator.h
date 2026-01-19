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
#include "includes/model_part.h"
#include "future/linear_operators/linear_operator.h"

namespace Kratos::Future
{

///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @brief
 * @class SparseMatrixLinearOperator
 * @ingroup KratosCore
 * @brief Linear operator handling a sparse matrix
 * @details
 * This class implements a LinearOperator that handles a sparse matrix,
 * allowing to use it in matrix-free algorithms while still storing
 * the matrix entries explicitly. The linear operator can be also used
 * in matrix-based algorithms by accessing the underlying sparse matrix.
 * @tparam TLinearAlgebra The struct containing the linear algebra types
 */
template <class TLinearAlgebra>
class SparseMatrixLinearOperator final : public LinearOperator<TLinearAlgebra>
{

public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SparseMatrixLinearOperator
    KRATOS_CLASS_POINTER_DEFINITION(SparseMatrixLinearOperator);

    /// Vector type definition from template parameter
    using VectorType = typename TLinearAlgebra::VectorType;

    /// Matrix type definition from template parameter
    using MatrixType = typename TLinearAlgebra::MatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @details Creates an empty SparseMatrixLinearOperator with uninitialized function objects.
     */
    SparseMatrixLinearOperator() = default;

    /**
     * @brief Constructor from a CSR matrix.
     * Constructs a SparseMatrixLinearOperator that wraps the provided CSR matrix.
     * Note that this constructor is only enabled when a CSR matrix type is provided.
     * @param rA Reference to the CSR matrix
     */
    SparseMatrixLinearOperator(MatrixType& rA)
        : LinearOperator<TLinearAlgebra>()
    {
        mrCsrMatrix = rA;
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
        const VectorType& rX,
        VectorType& rY) const override
    {
        mrCsrMatrix.SpMV(rX, rY);
    }

    void TransposeSpMV(
        const VectorType& rX,
        VectorType& rY) const override
    {
        mrCsrMatrix.TransposeSpMV(rX, rY);
    }

    ///@}
    ///@name Access
    ///@{

    MatrixType& GetMatrix()
    {
        return mrCsrMatrix;
    }

    const MatrixType& GetMatrix() const
    {
        return mrCsrMatrix;
    }

    ///@}
    ///@name Inquiry
    ///@{

    bool IsMatrixFree() const override
    {
        return false;
    }

    ///@}
private:

    /// Reference to the CSR matrix (of type MatrixType& handled as std::any)
    MatrixType& mrCsrMatrix;

}; // class SparseMatrixLinearOperator

///@}
///@name Input and output
///@{

///@}
///@} addtogroup block

} // namespace Kratos::Future
