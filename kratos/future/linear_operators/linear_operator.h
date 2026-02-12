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

namespace Kratos::Future
{

///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @brief
 * @class LinearOperator
 * @ingroup KratosCore
 * @brief Auxiliary container to store a linear operator
 * @details
 * This class implements a generic linear operator interface to be used
 * in matrix-free algorithms. Being a generic interface, it can be derived
 * to wrap different types of linear operators, including maxtrix-based ones.
 * @tparam TLinearAlgebra The struct containing the linear algebra types
 */
template <class TLinearAlgebra>
class LinearOperator
{

public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of LinearOperator
    KRATOS_CLASS_POINTER_DEFINITION(LinearOperator);

    /// Vector type definition from template parameter
    using VectorType = typename TLinearAlgebra::VectorType;

    /// Matrix type definition from template parameter
    using MatrixType = typename TLinearAlgebra::MatrixType;

    /// Data type stored in the system vector
    using DataType = typename VectorType::DataType;

    /// Index type used in the system vector
    using IndexType = typename VectorType::IndexType;

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
        mSize1 = ThisParameters["size_1"].GetInt();
        mSize2 = ThisParameters["size_2"].GetInt();
    }

    /**
     * @brief Construct a new Linear Operator object
     * @param Shape Pair with the shape of the linear operator
     */
    LinearOperator(const std::pair<std::size_t, std::size_t> Shape)
    {
        mSize1 = std::get<0>(Shape);
        mSize2 = std::get<1>(Shape);
    }

    /// Deleted copy constructor (non-copyable)
    LinearOperator(const LinearOperator& rOther) = delete;

    /// Defaulted move constructor
    LinearOperator(LinearOperator&& rOther) = default;

    /// Default destructor
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

    /**
     * @brief Performs the matrix-vector product y = A * x.
     * @param rX Input vector x
     * @param rY Output vector y
     */
    virtual void SpMV(
        const VectorType& rX,
        VectorType& rY) const
    {
        KRATOS_ERROR << "SpMV() is not implemented in base LinearOperator class." << std::endl;
    }

    /**
     * @brief Performs the transposed matrix-vector product y = A^T * x.
     * @param rX Input vector x
     * @param rY Output vector y
     */
    virtual void TransposeSpMV(
        const VectorType& rX,
        VectorType& rY) const
    {
        KRATOS_ERROR << "TransposeSpMV() is not implemented in base LinearOperator class." << std::endl;
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Get the Matrix object
     * @return MatrixType& Reference to the underlying matrix.
     */
    virtual MatrixType& GetMatrix()
    {
        KRATOS_ERROR << "GetMatrix() is not implemented in base LinearOperator class." << std::endl;
    }

    /**
     * @brief Get the Matrix object
     * @return const MatrixType& Reference to the underlying matrix.
     */
    virtual const MatrixType& GetMatrix() const
    {
        KRATOS_ERROR << "GetMatrix() is not implemented in base LinearOperator class." << std::endl;
    }

    ///@}
    ///@name Inquiry
    ///@{

    std::pair<std::size_t, std::size_t> Shape() const
    {
        return {mSize1, mSize2};
    }

    /**
     * @brief Get the number of rows.
     * @return Number of rows of the operator
     */
    std::size_t Size1() const
    {
        return mSize1;
    }

    /**
     * @brief Get the number of columns.
     * @return Number of columns of the operator
     */
    std::size_t Size2() const
    {
        return mSize2;
    }

    /**
     * @brief Check if the operator is matrix-free.
     * @return True if the operator is matrix-free, false otherwise
     */
    virtual bool IsMatrixFree() const
    {
        return true;
    }

    ///@}
private:

    ///@name Member Variables
    ///@{

    /// Number of rows of the operator
    std::size_t mSize1 = 0;

    /// Number of columns of the operator
    std::size_t mSize2 = 0;

    ///@}
}; // class LinearOperator

///@}
///@name Input and output
///@{

///@}
///@} addtogroup block

} // namespace Kratos::Future
