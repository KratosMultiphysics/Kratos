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
#include "future/linear_operators/sparse_matrix_linear_operator.h"

namespace Kratos::Future
{

///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

template <class TLinearAlgebra>
class LinearSystem final
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of LinearSystem
    KRATOS_CLASS_POINTER_DEFINITION(LinearSystem);

    /// DofsArrayType definition
    using DofsArrayType = ModelPart::DofsArrayType;

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
        std::string SystemName = "")
        : mpLhs(pLhs)
        , mpRhs(pRhs)
        , mpSol(pSol)
        , mSystemName(SystemName)
    {
        mpLinearOperator = Kratos::make_shared<SparseMatrixLinearOperator<TLinearAlgebra>>(*pLhs);
    }

    /// Constructor with matrix and additional data
    LinearSystem(
        typename MatrixType::Pointer pLhs,
        typename VectorType::Pointer pRhs,
        typename VectorType::Pointer pSol,
        ModelPart& rModelPart,
        DofsArrayType& rDofs,
        const std::string SystemName = "")
        : mpLhs(pLhs)
        , mpRhs(pRhs)
        , mpSol(pSol)
        , mSystemName(SystemName)
    {
        mpModelPart = &rModelPart;
        mpDofs = Kratos::make_shared<DofsArrayType>(rDofs);
        mpLinearOperator = Kratos::make_shared<SparseMatrixLinearOperator<TLinearAlgebra>>(*pLhs);
    }

    /// Constructor with linear operator
    LinearSystem(
        typename LinearOperatorType::Pointer pLinearOperator,
        typename VectorType::Pointer pRhs,
        typename VectorType::Pointer pSol,
        std::string SystemName = "")
        : mpLinearOperator(pLinearOperator)
        , mpRhs(pRhs)
        , mpSol(pSol)
        , mSystemName(SystemName)
    {
    }

    /// Constructor with linear operator and additional data
    LinearSystem(
        typename LinearOperatorType::Pointer pLinearOperator,
        typename VectorType::Pointer pRhs,
        typename VectorType::Pointer pSol,
        ModelPart& rModelPart,
        DofsArrayType& rDofs,
        const std::string SystemName = "")
        : mpLinearOperator(pLinearOperator)
        , mpRhs(pRhs)
        , mpSol(pSol)
        , mSystemName(SystemName)
    {
        mpModelPart = &rModelPart;
        mpDofs = Kratos::make_shared<DofsArrayType>(rDofs);
    }

    /// Copy constructor
    LinearSystem(const LinearSystem& rOther) = delete;

    /// Move constructor
    LinearSystem(LinearSystem&& rOther) = delete;

    /// Destructor
    ~LinearSystem() = default;

    ///@}
    ///@name Operators
    ///@{

    /// Copy assignment operator
    LinearSystem& operator=(const LinearSystem& rOther) = delete;

    /// Move assignment operator
    LinearSystem& operator=(LinearSystem&& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    virtual int Check()
    {
        // Check if the system arrays are initialized
        KRATOS_ERROR_IF(mpLinearOperator == nullptr && mpLhs == nullptr) << "Linear operator or left-hand side matrix must be initialized." << std::endl;
        KRATOS_ERROR_IF(mpRhs == nullptr) << "Right-hand side vector must be initialized." << std::endl;
        KRATOS_ERROR_IF(mpSol == nullptr) << "Solution vector must be initialized." << std::endl;

        // Check if the system name is set (not mandatory but a good practice for debugging purposes)
        KRATOS_WARNING_IF("LinearSystem", mSystemName.empty()) << "System name is empty." << std::endl;

        return 0;
    }

    ///@}
    ///@name Access
    ///@{

    /**
    * @brief Get a pointer to the linear operator.
    * @return Pointer to the linear operator
    */
    virtual typename LinearOperatorType::Pointer pGetLinearOperator()
    {
        KRATOS_ERROR_IF(!mpLinearOperator) << "Linear operator is not initialized." << std::endl;
        return mpLinearOperator;
    }

    /**
    * @brief Get a reference to the linear operator.
    * @return Reference to the linear operator
    */
    virtual LinearOperatorType& GetLinearOperator()
    {
        KRATOS_ERROR_IF(!mpLinearOperator) << "Linear operator is not initialized." << std::endl;
        return *mpLinearOperator;
    }

    /**
    * @brief Get a const reference to the linear operator.
    * @return Const reference to the linear operator
    */
    virtual const LinearOperatorType& GetLinearOperator() const
    {
        KRATOS_ERROR_IF(!mpLinearOperator) << "Linear operator is not initialized." << std::endl;
        return *mpLinearOperator;
    }

    /**
    * @brief Get a reference to the left-hand side matrix.
    * @return Reference to the left-hand side matrix
    */
    virtual typename MatrixType::Pointer pGetLeftHandSide()
    {
        KRATOS_ERROR_IF(!mpLhs) << "Left-hand side matrix is not initialized." << std::endl;
        return mpLhs;
    }

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
     * @brief Set the Left Hand Side matrix
     * @param rLhs The left-hand side matrix
     */
    virtual void SetLeftHandSide(MatrixType& rLhs)
    {
        mpLhs = Kratos::make_shared<MatrixType>(rLhs);
        mpLinearOperator = Kratos::make_shared<SparseMatrixLinearOperator<TLinearAlgebra>>(rLhs);
    }

    /**
    * @brief Get a pointer to the right-hand side vector.
    * @return Pointer to the right-hand side vector
    */
    virtual typename VectorType::Pointer pGetRightHandSide()
    {
        KRATOS_ERROR_IF(!mpRhs) << "Right-hand side vector is not initialized." << std::endl;
        return mpRhs;
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
     * @brief Set the Right Hand Side vector
     * @param rRhs The right-hand side vector
     */
    virtual void SetRightHandSide(VectorType& rRhs)
    {
        mpRhs = Kratos::make_shared<VectorType>(rRhs);
    }

    /**
     * @brief Get a pointer to the solution vector.
     * @return Pointer to the solution vector
     */
    virtual typename VectorType::Pointer pGetSolution() //FIXME: rename to pGetSolutionIncrement?
    {
        KRATOS_ERROR_IF(!mpSol) << "Solution vector is not initialized." << std::endl;
        return mpSol;
    }

    /**
    * @brief Get a reference to the solution vector.
    * @return Reference to the solution vector
    */
    virtual VectorType& GetSolution() //FIXME: rename to GetSolutionIncrement?
    {
        KRATOS_ERROR_IF(!mpSol) << "Solution vector is not initialized." << std::endl;
        return *mpSol;
    }

    /**
    * @brief Get a const reference to the solution vector.
    * @return Const reference to the solution vector
    */
    virtual const VectorType& GetSolution() const //FIXME: rename to GetSolutionIncrement?
    {
        KRATOS_ERROR_IF(!mpSol) << "Solution vector is not initialized." << std::endl;
        return *mpSol;
    }

    /**
     * @brief Set the Solution vector
     * @param rSol The solution vector
     */
    virtual void SetSolution(VectorType& rSol) //FIXME: rename to SetSolutionIncrement?
    {
        mpSol = Kratos::make_shared<VectorType>(rSol);
    }

    /**
     * @brief Set the Additional Data
     * This method is used to set the additional data of the linear system.
     * @param rModelPart The model part from which the linear system is built
     * @param rDofs The dofs array of the linear system
     */
    virtual void SetAdditionalData(
        ModelPart& rModelPart,
        DofsArrayType& rDofs)
    {
        KRATOS_WARNING_IF("LinearSystem", HasAdditionalData()) << "Additional data is already set for linear system " << Name() << ". Overwriting it." << std::endl;
        mpModelPart = &rModelPart;
        mpDofs = Kratos::make_shared<DofsArrayType>(rDofs);
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

    /**
     * @brief This method checks if the dimensions of the system of equations are consistent
     * @return True if consistent, false otherwise
     */
    virtual bool IsConsistent()
    {
        const std::size_t num_rows = mpLinearOperator->NumRows();
        const std::size_t num_cols = mpLinearOperator->NumCols();
        const std::size_t size_x = mpSol->size();
        const std::size_t size_b = mpRhs->size();
        return ((num_rows ==  num_cols) && (num_rows ==  size_x) && (num_rows ==  size_b));
    }

    /**
     * @brief This method checks if the dimensions of the system of equations are not consistent
     * @return False if consistent, true otherwise
     */
    virtual bool IsNotConsistent()
    {
        return !IsConsistent();
    }

    /**
    * @brief Get the name of the linear system.
    * @return Name of the linear system
    */
    virtual const std::string& Name() const
    {
        return mSystemName;
    }

    /**
     * @brief Check if the linear system has additional data (model part and dofs)
     * @return true if the linear system has additional data
     * @return false if the linear system does not have additional data
     */
    virtual bool HasAdditionalData() const
    {
        return mpModelPart != nullptr && mpDofs != nullptr;
    }

    ///@}
private:

    ///@name Member Variables
    ///@{

    typename LinearOperatorType::Pointer mpLinearOperator = nullptr;

    typename MatrixType::Pointer mpLhs = nullptr;

    typename VectorType::Pointer mpRhs = nullptr;

    typename VectorType::Pointer mpSol = nullptr;

    const std::string mSystemName;

    ModelPart* mpModelPart = nullptr;

    typename ModelPart::DofsArrayType::Pointer mpDofs = nullptr;

    ///@}
}; // Class LinearSystem

///@}
///@name Input and output
///@{

///@}
///@} addtogroup block

} // namespace Kratos.