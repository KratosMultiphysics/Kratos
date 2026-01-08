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

    /**
     * @brief Performs the matrix-vector product y = A * x.
     * @param rX Input vector x
     * @param rY Output vector y
     */
    virtual void SpMV(
        const VectorType& rX,
        VectorType& rY) const = 0;

    /**
     * @brief Performs the transposed matrix-vector product y = A^T * x.
     * @param rX Input vector x
     * @param rY Output vector y
     */
    virtual void TransposeSpMV(
        const VectorType& rX,
        VectorType& rY) const = 0;

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

    /**
     * @brief Get a reference to the underlying matrix.
     * @tparam TMatrixType The type of the matrix to retrieve
     * @return Reference to the matrix
     */
    template<class TMatrixType>
    TMatrixType& GetMatrix()
    {
        auto& r_matrix = this->GetMatrixImpl(); // Get the underlying matrix as std::any
        return std::any_cast<std::reference_wrapper<TMatrixType>>(r_matrix).get(); // Cast and return the reference
    }

    /**
     * @brief Get a const reference to the underlying matrix.
     * @tparam TMatrixType The type of the matrix to retrieve
     * @return Const reference to the matrix
     */
    template<class TMatrixType>
    const TMatrixType& GetMatrix() const
    {
        const auto& r_matrix = this->GetMatrixImpl(); // Get the underlying matrix as const std::any
        return std::any_cast<std::reference_wrapper<TMatrixType>>(r_matrix).get(); // Cast and return the reference
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
     * @brief Implementation to get a reference to the underlying matrix.
     * This method is to be overridden in matrix-based LinearOperator classes.
     * An exception is thrown if this method is called from the base class (matrix-free).
     * @return Reference to the matrix as an std::any
     */
    [[noreturn]] virtual std::any& GetMatrixImpl()
    {
        KRATOS_ERROR << "GetMatrixImpl() not implemented in base LinearOperator class." << std::endl;
    }

    /**
     * @brief Implementation to get a const reference to the underlying matrix.
     * This method is to be overridden in matrix-based LinearOperator classes.
     * An exception is thrown if this method is called from the base class (matrix-free).
     * @return Reference to the matrix as a const std::any
     */
    [[noreturn]] virtual const std::any& GetMatrixImpl() const
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
