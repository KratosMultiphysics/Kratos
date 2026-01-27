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
#include "future/linear_operators/linear_operator.h"

namespace Kratos::Future
{

///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

template <class TLinearAlgebra>
class EigenvalueSystem final
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

    /// Type od the dense matrix
    using DenseMatrixType = typename TLinearAlgebra::DenseMatrixType;

    /// Type of the dense matrix pointer
    using DenseMatrixPointerType = Kratos::shared_ptr<DenseMatrixType>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    EigenvalueSystem()
    {
    }

    /// Constructor with matrix
    EigenvalueSystem(
        typename MatrixType::Pointer pK,
        typename MatrixType::Pointer pM,
        typename VectorType::Pointer pEigenvalues,
        DenseMatrixPointerType pEigenvectors,
        const std::string SystemName = "")
        : mpK(pK)
        , mpM(pM)
        , mpEigenvalues(pEigenvalues)
        , mpEigenvectors(pEigenvectors)
        , mSystemName(SystemName)
    {
    }

    /// Copy constructor
    EigenvalueSystem(const EigenvalueSystem& rOther) = delete;

    /// Move constructor
    EigenvalueSystem(EigenvalueSystem&& rOther) = delete;

    /// Destructor
    ~EigenvalueSystem() = default;

    ///@}
    ///@name Operators
    ///@{

    /// Copy assignment operator
    EigenvalueSystem& operator=(const EigenvalueSystem& rOther) = delete;

    /// Defaulted move assignment operator
    EigenvalueSystem& operator=(EigenvalueSystem&& rOther) = delete;

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
    DenseMatrixType& GetEigenvectors()
    {
        KRATOS_ERROR_IF(!mpEigenvectors) << "Eigenvectors matrix is not initialized." << std::endl;
        return *mpEigenvectors;
    }

    /**
    * @brief Get a const reference to the eigenvectors matrix.
    * @return Const reference to the eigenvectors matrix
    */
    const DenseMatrixType& GetEigenvectors() const
    {
        KRATOS_ERROR_IF(!mpEigenvectors) << "Eigenvectors matrix is not initialized." << std::endl;
        return *mpEigenvectors;
    }

    ///@}
    ///@name Inquiry
    ///@{

    /**
    * @brief Get the name of the linear system.
    * @return Name of the linear system
    */
    virtual const std::string& Name() const
    {
        return mSystemName;
    }

    ///@}
private:

    ///@name Member Variables
    ///@{

    typename MatrixType::Pointer mpK = nullptr;

    typename MatrixType::Pointer mpM = nullptr;

    typename VectorType::Pointer mpEigenvalues = nullptr;

    DenseMatrixPointerType mpEigenvectors = nullptr;

    const std::string mSystemName;

    ///@}
}; // Class EigenvalueSystem

///@}
///@name Input and output
///@{

///@}
///@} addtogroup block

} // namespace Kratos.