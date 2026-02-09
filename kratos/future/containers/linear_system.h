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

/// Enum defining the linar system arrays tags
enum class LinearSystemTag
{
    LHS = 0, // Left hand side matrix tag
    RHS = 1, // Right hand side vector tag
    Dx = 2, // Solution increment tag
    M = 3, // Mass matrix tag
    K = 4, // Stiffness matrix tag
    C = 5, // Damping matrix tag
    Eigvals = 6, // Eigenvalues tag
    Eigvects = 7, // Eigenvectors tag
    NumberOfTags = 8 // Sentinel tag with the size
};

/**
 * @brief Linear system container
 * This class encapsulates the components of a linear system, including the
 * left-hand side matrix (either in matrix-based or matrix-free form), the
 * right-hand side vector and the solution vector.
 * The storage of the matrices is handled by a linear operator.
 * The ownership of the linear operator is transferred to the linear system but not that of the arrays.
 * Information about the physics of the problem can be optionally stored via a ModelPart and its Dofs.
 * @tparam TLinearAlgebra The struct containing the linear algebra types
 */
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

    /// Default onstructor
    LinearSystem() = default;

    /// Constructor with system name only
    LinearSystem(std::string SystemName)
        : mSystemName(SystemName)
    {
    }

    /// Constructor with left hand side matrix
    LinearSystem(
        typename MatrixType::Pointer pLhs,
        typename VectorType::Pointer pRhs,
        typename VectorType::Pointer pDx,
        std::string SystemName = "")
        : mSystemName(SystemName)
    {
        this->pSetMatrix(pLhs, LinearSystemTag::LHS);
        this->pSetVector(pRhs, LinearSystemTag::RHS);
        this->pSetVector(pDx, LinearSystemTag::Dx);
    }

    /// Constructor with linear operator emulating the left hand side matrix
    LinearSystem(
        typename LinearOperatorType::UniquePointer&& pLinearOperator,
        typename VectorType::Pointer pRhs,
        typename VectorType::Pointer pDx,
        std::string SystemName = "")
        : mSystemName(SystemName)
    {
        this->pSetLinearOperator(std::move(pLinearOperator), LinearSystemTag::LHS);
        this->pSetVector(pRhs, LinearSystemTag::RHS);
        this->pSetVector(pDx, LinearSystemTag::Dx);
    }

    /// Constructor with matrix and additional data
    LinearSystem(
        typename MatrixType::Pointer pLhs,
        typename VectorType::Pointer pRhs,
        typename VectorType::Pointer pDx,
        ModelPart& rModelPart,
        DofsArrayType& rDofs,
        const std::string SystemName = "")
        : mSystemName(SystemName)
    {
        mpModelPart = &rModelPart;
        mpDofs = Kratos::make_shared<DofsArrayType>(rDofs);
        this->pSetMatrix(pLhs, LinearSystemTag::LHS);
        this->pSetVector(pRhs, LinearSystemTag::RHS);
        this->pSetVector(pDx, LinearSystemTag::Dx);
    }

    /// Constructor with linear operator and additional data
    LinearSystem(
        typename LinearOperatorType::Pointer pLinearOperator,
        typename VectorType::Pointer pRhs,
        typename VectorType::Pointer pSol,
        ModelPart& rModelPart,
        DofsArrayType& rDofs,
        const std::string SystemName = "")
        : mSystemName(SystemName)
    {
        mpModelPart = &rModelPart;
        mpDofs = Kratos::make_shared<DofsArrayType>(rDofs);
        this->pSetLinearOperator(std::move(pLinearOperator), LinearSystemTag::LHS);
        this->pSetVector(pRhs, LinearSystemTag::RHS);
        this->pSetVector(pSol, LinearSystemTag::Dx);
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

    /**
     * @brief This method checks if the linear system is consistent
     * @return 0 if consistent
     */
    virtual int Check()
    {
        // Check if the system name is set (not mandatory but a good practice for debugging purposes)
        KRATOS_WARNING_IF("LinearSystem", mSystemName.empty()) << "System name is empty." << std::endl;

        return 0;
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Get a matrix
     * @param Tag The tag of the matrix to be retrieved (@see LinearSystemTag)
     * @return Reference to the matrix
     */
    MatrixType& GetMatrix(LinearSystemTag Tag) const
    {
        const std::size_t tag_idx = GetTagIndex(Tag);
        const auto& p_matrix_lin_op = mLinearOperators[tag_idx];
        KRATOS_ERROR_IF(!p_matrix_lin_op) << "Requested matrix is not initialized for tag " << tag_idx << std::endl;
        KRATOS_ERROR_IF(p_matrix_lin_op->IsMatrixFree()) << "Linear operator referring to matrix " << tag_idx << " is matrix free" << std::endl;
        return p_matrix_lin_op->GetMatrix();
    }

    /**
     * @brief Set a matrix to an specific tag
     * @details Note that the matrices are stored by creating a matrix-based linear operator
     * @param pMatrix Pointer to the matrix
     * @param Tag The tag of the matrix to be set (@see LinearSystemTag)
     */
    void pSetMatrix(
        typename MatrixType::Pointer pMatrix,
        LinearSystemTag Tag)
    {
        KRATOS_ERROR_IF(!pMatrix) << "Provided matrix pointer is null" << std::endl;
        mLinearOperators[GetTagIndex(Tag)] = Kratos::make_unique<SparseMatrixLinearOperator<TLinearAlgebra>>(*pMatrix);
    }

    /**
     * @brief Get a vector
     * @param Tag The tag of the vector to be retrieved (@see LinearSystemTag)
     * @return Reference to the vector
     */
    VectorType& GetVector(LinearSystemTag Tag) const
    {
        const std::size_t tag_idx = GetTagIndex(Tag);
        const auto& p_vector = mVectors[tag_idx];
        KRATOS_ERROR_IF(!p_vector) << "Requested vector is not initialized for tag " << tag_idx << std::endl;
        return *p_vector;
    }

    /**
     * @brief Set a vector to an specific tag
     * @param pVector Pointer to the vector
     * @param Tag The tag of the vector to be set (@see LinearSystemTag)
     */
    void pSetVector(
        typename VectorType::Pointer pVector,
        LinearSystemTag Tag)
    {
        KRATOS_ERROR_IF(!pVector) << "Provided vector pointer is null" << std::endl;
        mVectors[GetTagIndex(Tag)] = pVector;
    }

    /**
     * @brief Get a linear operator
     * @param Tag The tag of the linear operator to be retrieved (@see LinearSystemTag)
     * @return Reference to the linear operator
     */
    LinearOperatorType& GetLinearOperator(LinearSystemTag Tag) const
    {
        const std::size_t tag_idx = GetTagIndex(Tag);
        const auto& p_linear_operator = mLinearOperators[tag_idx];
        KRATOS_ERROR_IF(!p_linear_operator) << "Requested linear operator is not initialized for tag " << tag_idx << std::endl;
        return *p_linear_operator;
    }

    /**
     * @brief Set a linear operator to an specific tag
     * @details Note that the linear operator ownership is transferred to the linear system
     * so a linear operator cannot be set to more than one linear system or tag.
     * @param pLinearOperator Pointer to the linear operator
     * @param Tag The tag of the linear operator to be set (@see LinearSystemTag)
     */
    void pSetLinearOperator(
        typename LinearOperatorType::UniquePointer&& pLinearOperator,
        LinearSystemTag Tag)
    {
        KRATOS_ERROR_IF(!pLinearOperator) << "Provided linear operator pointer is null" << std::endl;
        mLinearOperators[GetTagIndex(Tag)] = std::move(pLinearOperator);
    }

    /**
     * @brief Set the Additional Data
     * This method is used to set the additional data of the linear system.
     * @param rModelPart The model part from which the linear system is built
     * @param rDofs The dofs array of the linear system
     */
    void SetAdditionalData(
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
     * @brief This method checks if the dimensions of the system of equations are consistent
     * @return True if consistent, false otherwise
     */
    bool IsConsistent()
    {
        const std::size_t num_rows = GetLinearOperator(LinearSystemTag::LHS).NumRows();
        const std::size_t num_cols = GetLinearOperator(LinearSystemTag::LHS).NumCols();
        const std::size_t size_x = GetVector(LinearSystemTag::Dx).size();
        const std::size_t size_b = GetVector(LinearSystemTag::RHS).size();
        return ((num_rows ==  num_cols) && (num_rows ==  size_x) && (num_rows ==  size_b));
    }

    /**
     * @brief This method checks if the dimensions of the system of equations are not consistent
     * @return False if consistent, true otherwise
     */
    bool IsNotConsistent()
    {
        return !IsConsistent();
    }

    /**
    * @brief Get the name of the linear system.
    * @return Name of the linear system
    */
    const std::string& Name() const
    {
        return mSystemName;
    }

    /**
     * @brief Check if the linear system has additional data (model part and dofs)
     * @return true if the linear system has additional data
     * @return false if the linear system does not have additional data
     */
    bool HasAdditionalData() const
    {
        return mpModelPart != nullptr && mpDofs != nullptr;
    }

    ///@}
private:

    ///@name Member Variables
    ///@{

    const std::string mSystemName; // Name of the linear system

    ModelPart* mpModelPart = nullptr; // Model part of the linear system

    typename ModelPart::DofsArrayType::Pointer mpDofs = nullptr; // Dofs of the linear system

    std::array<typename VectorType::Pointer, static_cast<std::size_t>(LinearSystemTag::NumberOfTags)> mVectors; // Vectors of the linear system

    std::array<typename LinearOperatorType::UniquePointer, static_cast<std::size_t>(LinearSystemTag::NumberOfTags)> mLinearOperators; // Linear operators of the linear system

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief Get the index of a tag
     * @param Tag The tag of the linear operator to be retrieved (@see LinearSystemTag)
     * @return constexpr std::size_t The index of the tag
     */
    static constexpr std::size_t GetTagIndex(LinearSystemTag Tag)
    {
        return static_cast<std::size_t>(Tag);
    }

    ///@}
}; // Class LinearSystem

///@}
///@name Input and output
///@{

///@}
///@} addtogroup block

} // namespace Kratos.