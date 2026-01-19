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
#include "future/containers/linear_system.h"
#include "future/linear_operators/linear_operator.h"

namespace Kratos::Future
{

///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

template <class TLinearAlgebra>
class EigenvalueSystem : public LinearSystem<TLinearAlgebra>
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of EigenvalueSystem
    KRATOS_CLASS_POINTER_DEFINITION(EigenvalueSystem);

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
    EigenvalueSystem()
        : LinearSystem<TLinearAlgebra>()
    {
    }

    /// Constructor with matrix
    EigenvalueSystem(
        typename MatrixType::Pointer pK,
        typename MatrixType::Pointer pM,
        typename VectorType::Pointer pEigenvalues,
        typename VectorType::Pointer pEigenvectors)
        : LinearSystem<TLinearAlgebra>()
        , mpK(pK)
        , mpM(pM)
        , mpEigenvalues(pEigenvalues)
        , mpEigenvectors(pEigenvectors)
    {
    }

    // /// Copy constructor
    // EigenvalueSystem(const EigenvalueSystem& rOther)
    // {
    //     mpLhs = rOther.mpLhs;
    //     mpRhs = rOther.mpRhs;
    //     mpLinearOperator = rOther.mpLinearOperator;
    // }

    // /// Defaulted move constructor
    // EigenvalueSystem(EigenvalueSystem&& rOther)
    // {
    //     mpLhs = std::move(rOther.mpLhs);
    //     mpRhs = std::move(rOther.mpRhs);
    //     mpLinearOperator = std::move(rOther.mpLinearOperator);
    // }

    /// Destructor
    ~EigenvalueSystem() = default;

    ///@}
    ///@name Operators
    ///@{

    // /// Copy assignment operator
    // EigenvalueSystem& operator=(const EigenvalueSystem& rOther)
    // {
    //     mpLhs = rOther.mpLhs;
    //     mpRhs = rOther.mpRhs;
    //     mpLinearOperator = rOther.mpLinearOperator;
    //     return *this;
    // }

    // /// Defaulted move assignment operator
    // EigenvalueSystem& operator=(EigenvalueSystem&& rOther)
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
    * @brief Get a reference to the stiffness matrix.
    * @return Reference to the stiffness matrix
    */
    MatrixType& GetStiffnessMatrix()
    {
        KRATOS_ERROR_IF(!mpK) << "Stiffness matrix is not initialized." << std::endl;
        return *mpK;
    }

    /**
    * @brief Get a const reference to the stiffness matrix.
    * @return Const reference to the stiffness matrix
    */
    const MatrixType& GetStiffnessMatrix() const
    {
        KRATOS_ERROR_IF(!mpK) << "Stiffness matrix is not initialized." << std::endl;
        return *mpK;
    }

    /**
    * @brief Get a reference to the mass matrix.
    * @return Reference to the mass matrix
    */
    MatrixType& GetMassMatrix()
    {
        KRATOS_ERROR_IF(!mpM) << "Mass matrix is not initialized." << std::endl;
        return *mpM;
    }

    /**
    * @brief Get a const reference to the mass matrix.
    * @return Const reference to the mass matrix
    */
    const MatrixType& GetMassMatrix() const
    {
        KRATOS_ERROR_IF(!mpM) << "Mass matrix is not initialized." << std::endl;
        return *mpM;
    }

    /**
    * @brief Get a reference to the eigenvalues vector.
    * @return Reference to the eigenvalues vector
    */
    VectorType& GetEigenvalues()
    {
        KRATOS_ERROR_IF(!mpEigenvalues) << "Eigenvalues vector is not initialized." << std::endl;
        return *mpEigenvalues;
    }

    /**
    * @brief Get a const reference to the eigenvalues vector.
    * @return Const reference to the eigenvalues vector
    */
    const VectorType& GetEigenvalues() const
    {
        KRATOS_ERROR_IF(!mpEigenvalues) << "Eigenvalues vector is not initialized." << std::endl;
        return *mpEigenvalues;
    }

    /**
    * @brief Get a reference to the eigenvectors matrix.
    * @return Reference to the eigenvectors matrix
    */
    MatrixType& GetEigenvectors()
    {
        KRATOS_ERROR_IF(!mpEigenvectors) << "Eigenvectors matrix is not initialized." << std::endl;
        return *mpEigenvectors;
    }

    /**
    * @brief Get a const reference to the eigenvectors matrix.
    * @return Const reference to the eigenvectors matrix
    */
    const MatrixType& GetEigenvectors() const
    {
        KRATOS_ERROR_IF(!mpEigenvectors) << "Eigenvectors matrix is not initialized." << std::endl;
        return *mpEigenvectors;
    }

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
private:

    ///@name Member Variables
    ///@{

    typename MatrixType::Pointer mpK = nullptr;

    typename MatrixType::Pointer mpM = nullptr;

    typename VectorType::Pointer mpEigenvalues = nullptr;

    typename MultipleVectorType::Pointer mpEigenvectors = nullptr;

    ///@}
}; // Class EigenvalueSystem

///@}
///@name Input and output
///@{

///@}
///@} addtogroup block

} // namespace Kratos.