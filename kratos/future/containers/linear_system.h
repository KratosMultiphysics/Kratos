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
//                   Riccardo Rossi
//

#pragma once

// System includes

// External includes

// Project includes
#include "future/linear_operators/linear_operator.h"
#include "future/linear_operators/sparse_matrix_linear_operator.h"

namespace Kratos::Future
{

///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

template <class TLinearAlgebra>
class LinearSystem
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of LinearSystem
    KRATOS_CLASS_POINTER_DEFINITION(LinearSystem);

    /// Type of the matrix from the linear algebra traits
    using MatrixType = typename TLinearAlgebra::MatrixType;

    /// Type of the vector from the linear algebra traits
    using VectorType = typename TLinearAlgebra::VectorType;

    /// Type of the linear operator from the linear algebra traits
    using LinearOperatorType = LinearOperator<TLinearAlgebra>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    LinearSystem() = default;

    /// Constructor with matrix
    LinearSystem(
        typename MatrixType::Pointer pLhs,
        typename VectorType::Pointer pRhs,
        typename VectorType::Pointer pSol,
        std::string system_name = "")
        : mpLhs(pLhs)
        , mpRhs(pRhs)
        , mpSol(pSol)
    {
        mpLinearOperator = Kratos::make_shared<SparseMatrixLinearOperator<TLinearAlgebra>>(*pLhs);
    }

    /// Constructor with linear operator
    LinearSystem(
        typename LinearOperatorType::Pointer pLinearOperator,
        typename VectorType::Pointer pRhs,
        typename VectorType::Pointer pSol)
        : mpLinearOperator(pLinearOperator)
        , mpRhs(pRhs)
        , mpSol(pSol)
    {
    }

    // /// Copy constructor
    // LinearSystem(const LinearSystem& rOther)
    // {
    //     mpLhs = rOther.mpLhs;
    //     mpRhs = rOther.mpRhs;
    //     mpLinearOperator = rOther.mpLinearOperator;
    // }

    // /// Defaulted move constructor
    // LinearSystem(LinearSystem&& rOther)
    // {
    //     mpLhs = std::move(rOther.mpLhs);
    //     mpRhs = std::move(rOther.mpRhs);
    //     mpLinearOperator = std::move(rOther.mpLinearOperator);
    // }

    /// Destructor
    ~LinearSystem() = default;

    ///@}
    ///@name Operators
    ///@{

    // /// Copy assignment operator
    // LinearSystem& operator=(const LinearSystem& rOther)
    // {
    //     mpLhs = rOther.mpLhs;
    //     mpRhs = rOther.mpRhs;
    //     mpLinearOperator = rOther.mpLinearOperator;
    //     return *this;
    // }

    // /// Defaulted move assignment operator
    // LinearSystem& operator=(LinearSystem&& rOther)
    // {
    //     mpLhs = std::move(rOther.mpLhs);
    //     mpRhs = std::move(rOther.mpRhs);
    //     mpLinearOperator = std::move(rOther.mpLinearOperator);
    //     return *this;
    // }

    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Access
    ///@{

    /**
    * @brief Get a reference to the left-hand side matrix.
    * @return Reference to the left-hand side matrix
    */
    virtual MatrixType& GetLeftHandSide()
    {
        KRATOS_ERROR_IF(!mpLhs) << "Left-hand side matrix is not initialized." << std::endl;
        return *mpLhs;
    }

    /**
    * @brief Get a const reference to the left-hand side matrix.
    * @return Const reference to the left-hand side matrix
    */
    virtual const MatrixType& GetLeftHandSide() const
    {
        KRATOS_ERROR_IF(!mpLhs) << "Left-hand side matrix is not initialized." << std::endl;
        return *mpLhs;
    }

    /**
    * @brief Get a reference to the right-hand side vector.
    * @return Reference to the right-hand side vector
    */
    virtual VectorType& GetRightHandSide()
    {
        KRATOS_ERROR_IF(!mpRhs) << "Right-hand side vector is not initialized." << std::endl;
        return *mpRhs;
    }

    /**
    * @brief Get a const reference to the right-hand side vector.
    * @return Const reference to the right-hand side vector
    */
    virtual const VectorType& GetRightHandSide() const
    {
        KRATOS_ERROR_IF(!mpRhs) << "Right-hand side vector is not initialized." << std::endl;
        return *mpRhs;
    }

    /**
    * @brief Get a reference to the solution vector.
    * @return Reference to the solution vector
    */
    virtual VectorType& GetSolution()
    {
        KRATOS_ERROR_IF(!mpSol) << "Solution vector is not initialized." << std::endl;
        return *mpSol;
    }

    /**
    * @brief Get a const reference to the solution vector.
    * @return Const reference to the solution vector
    */
    virtual const VectorType& GetSolution() const
    {
        KRATOS_ERROR_IF(!mpSol) << "Solution vector is not initialized." << std::endl;
        return *mpSol;
    }

    ///@}
    ///@name Inquiry
    ///@{

    /**
    * @brief Check if the linear system is matrix-free.
    * @return True if the linear system is matrix-free, false if it is not
    */
    virtual bool IsMatrixFree() const
    {
        if (mpLhs) {
            return false;
        } else if (mpLinearOperator) {
            return true;
        } else {
            KRATOS_ERROR << "Linear system has no matrix or linear operator." << std::endl;
        }
    }

    ///@}
private:

    ///@name Member Variables
    ///@{

    typename VectorType::Pointer mpSol = nullptr;

    typename VectorType::Pointer mpRhs = nullptr;

    typename MatrixType::Pointer mpLhs = nullptr;

    typename LinearOperatorType::Pointer mpLinearOperator = nullptr;

    ///@}
}; // Class LinearSystem

///@}
///@name Input and output
///@{

///@}
///@} addtogroup block

} // namespace Kratos.