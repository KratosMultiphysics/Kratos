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
#include "future/containers/linear_system_tags.h"
#include "future/linear_operators/linear_operator.h"
#include "future/linear_operators/sparse_matrix_linear_operator.h"

namespace Kratos::Future
{

///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

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

    /// Type of the dense matrix from the linear algebra traits
    using DenseMatrixType = typename TLinearAlgebra::DenseMatrixType;

    /// Type of the dense matrix pointer
    using DenseMatrixPointerType = Kratos::shared_ptr<DenseMatrixType>;

    /// Type of the linear operator from the linear algebra traits
    using LinearOperatorType = LinearOperator<TLinearAlgebra>;

    /// Dense vector tag type definition
    using DenseVectorTag = typename LinearSystemTags::DenseVectorTag;

    /// Dense matrix tag type definition
    using DenseMatrixTag = typename LinearSystemTags::DenseMatrixTag;

    /// Sparse matrix tag type definition
    using SparseMatrixTag = typename LinearSystemTags::SparseMatrixTag;

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
        this->pSetMatrix(pLhs, SparseMatrixTag::LHS);
        this->pSetVector(pRhs, DenseVectorTag::RHS);
        this->pSetVector(pDx, DenseVectorTag::Dx);
    }

    /// Constructor with linear operator emulating the left hand side matrix
    LinearSystem(
        typename LinearOperatorType::UniquePointer&& pLinearOperator,
        typename VectorType::Pointer pRhs,
        typename VectorType::Pointer pDx,
        std::string SystemName = "")
        : mSystemName(SystemName)
    {
        this->pSetLinearOperator(std::move(pLinearOperator), SparseMatrixTag::LHS);
        this->pSetVector(pRhs, DenseVectorTag::RHS);
        this->pSetVector(pDx, DenseVectorTag::Dx);
    }

    /// Constructor with matrix and additional data
    LinearSystem(
        typename MatrixType::Pointer pLhs,
        typename VectorType::Pointer pRhs,
        typename VectorType::Pointer pDx,
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const std::string SystemName = "")
        : mSystemName(SystemName)
    {
        mpDofSet = &rDofSet;
        mpModelPart = &rModelPart;
        this->pSetMatrix(pLhs, SparseMatrixTag::LHS);
        this->pSetVector(pRhs, DenseVectorTag::RHS);
        this->pSetVector(pDx, DenseVectorTag::Dx);
    }

    /// Constructor with linear operator and additional data
    LinearSystem(
        typename LinearOperatorType::Pointer pLinearOperator,
        typename VectorType::Pointer pRhs,
        typename VectorType::Pointer pSol,
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const std::string SystemName = "")
        : mSystemName(SystemName)
    {
        mpDofSet = &rDofSet;
        mpModelPart = &rModelPart;
        this->pSetLinearOperator(std::move(pLinearOperator), SparseMatrixTag::LHS);
        this->pSetVector(pRhs, DenseVectorTag::RHS);
        this->pSetVector(pSol, DenseVectorTag::Dx);
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
     * @param Tag The tag of the matrix to be retrieved (@see SparseMatrixTag)
     * @return Pointer to the matrix
     */
    typename MatrixType::Pointer pGetMatrix(SparseMatrixTag Tag) const
    {
        const std::size_t tag_idx = LinearSystemTags::GetTagIndex(Tag);
        const auto& p_matrix_lin_op = mLinearOperators[tag_idx];
        KRATOS_ERROR_IF(!p_matrix_lin_op) << "Requested matrix is not initialized for tag " << tag_idx << std::endl;
        KRATOS_ERROR_IF(p_matrix_lin_op->IsMatrixFree()) << "Linear operator referring to matrix " << tag_idx << " is matrix free" << std::endl;
        return p_matrix_lin_op->pGetMatrix();
    }

    /**
     * @brief Set a matrix to an specific tag
     * @details Note that the matrices are stored by creating a matrix-based linear operator
     * @param pMatrix Pointer to the matrix
     * @param Tag The tag of the matrix to be set (@see SparseMatrixTag)
     */
    void pSetMatrix(
        typename MatrixType::Pointer pMatrix,
        SparseMatrixTag Tag)
    {
        KRATOS_ERROR_IF(!pMatrix) << "Provided matrix pointer is null" << std::endl;
        mLinearOperators[LinearSystemTags::GetTagIndex(Tag)] = Kratos::make_unique<SparseMatrixLinearOperator<TLinearAlgebra>>(pMatrix);
    }

    /**
     * @brief Get a dense matrix
     * @param Tag The tag of the dense matrix to be retrieved (@see DenseMatrixTag)
     * @return Pointer to the dense matrix
     */
    DenseMatrixPointerType pGetMatrix(DenseMatrixTag Tag) const
    {
        const std::size_t tag_idx = LinearSystemTags::GetTagIndex(Tag);
        const auto& p_dense_matrix = mDenseMatrices[tag_idx];
        KRATOS_ERROR_IF(!p_dense_matrix) << "Requested dense matrix is not initialized for tag " << tag_idx << std::endl;
        return p_dense_matrix;
    }

    /**
     * @brief Set a dense matrix to an specific tag
     * @param pDenseMatrix Pointer to the dense matrix
     * @param Tag The tag of the dense matrix to be set (@see DenseMatrixTag)
     */
    void pSetMatrix(
        DenseMatrixPointerType pDenseMatrix,
        DenseMatrixTag Tag)
    {
        KRATOS_ERROR_IF(!pDenseMatrix) << "Provided dense matrix pointer is null" << std::endl;
        mDenseMatrices[LinearSystemTags::GetTagIndex(Tag)] = pDenseMatrix;
    }

    /**
     * @brief Get a vector
     * @param Tag The tag of the vector to be retrieved (@see DenseVectorTag)
     * @return Pointer to the vector
     */
    typename VectorType::Pointer pGetVector(DenseVectorTag Tag) const
    {
        const std::size_t tag_idx = LinearSystemTags::GetTagIndex(Tag);
        const auto& p_vector = mVectors[tag_idx];
        KRATOS_ERROR_IF(!p_vector) << "Requested vector is not initialized for tag " << tag_idx << std::endl;
        return p_vector;
    }

    /**
     * @brief Set a vector to an specific tag
     * @param pVector Pointer to the vector
     * @param Tag The tag of the vector to be set (@see DenseVectorTag)
     */
    void pSetVector(
        typename VectorType::Pointer pVector,
        DenseVectorTag Tag)
    {
        KRATOS_ERROR_IF(!pVector) << "Provided vector pointer is null" << std::endl;
        mVectors[LinearSystemTags::GetTagIndex(Tag)] = pVector;
    }

    /**
     * @brief Get a linear operator
     * @param Tag The tag of the linear operator to be retrieved (@see SparseMatrixTag)
     * @return Pointer to the linear operator
     */
    const typename LinearOperatorType::UniquePointer& pGetLinearOperator(SparseMatrixTag Tag) const
    {
        const std::size_t tag_idx = LinearSystemTags::GetTagIndex(Tag);
        const auto& p_linear_operator = mLinearOperators[tag_idx];
        KRATOS_ERROR_IF(!p_linear_operator) << "Requested linear operator is not initialized for tag " << tag_idx << std::endl;
        return p_linear_operator;
    }

    /**
     * @brief Set a linear operator to an specific tag
     * @details Note that the linear operator ownership is transferred to the linear system
     * so a linear operator cannot be set to more than one linear system or tag.
     * @param pLinearOperator Pointer to the linear operator
     * @param Tag The tag of the linear operator to be set (@see SparseMatrixTag)
     */
    void pSetLinearOperator(
        typename LinearOperatorType::UniquePointer&& pLinearOperator,
        SparseMatrixTag Tag)
    {
        KRATOS_ERROR_IF(!pLinearOperator) << "Provided linear operator pointer is null" << std::endl;
        mLinearOperators[LinearSystemTags::GetTagIndex(Tag)] = std::move(pLinearOperator);
    }

    /**
     * @brief Set the Additional Data
     * @details Some solvers may require a minimum degree of knowledge of the structure of the matrix.
     * For example, when solving a mixed u-p problem, it is important to identify the row associated with v and p.
     * Another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers,
     * which require knowledge of the spatial position of the nodes associated with a given degree of freedom (DOF).
     * This function provides the opportunity to provide such data if needed.
     * @param rModelPart The model part from which the linear system is built
     * @param rDofSet The dofs array of the linear system
     */
    void SetAdditionalData(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet)
    {
        KRATOS_WARNING_IF("LinearSystem", HasAdditionalData()) << "Additional data is already set for linear system " << Name() << ". Overwriting it" << std::endl;
        mpDofSet = &rDofSet;
        mpModelPart = &rModelPart;
    }

    /**
     * @brief Get the Model Part
     * @return Pointer to the model part
     */
    ModelPart* pGetModelPart() const
    {
        KRATOS_ERROR_IF_NOT(mpModelPart) << "Model part is not set for linear system " << Name() << std::endl;
        return mpModelPart;
    }

    /**
     * @brief Get the Dofs Array
     * @return Pointer to the dofs array
     */
    DofsArrayType* pGetDofSet() const
    {
        KRATOS_ERROR_IF_NOT(mpDofSet) << "DOFs array is not set for linear system " << Name() << std::endl;
        return mpDofSet;
    }

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief This method checks if the dimensions of the given linear system are consistent
     * @param LhsMatrixTag The tag of the system matrix to be checked (@see SparseMatrixTag)
     * @param RhsVectorTag The tag of the residual vector to be checked (@see DenseVectorTag)
     * @param DxVectorTag The tag of the solution vector to be checked (@see DenseVectorTag)
     * @return True if consistent, false otherwise
     */
    bool IsConsistent(
        SparseMatrixTag LhsMatrixTag,
        DenseVectorTag RhsVectorTag,
        DenseVectorTag DxVectorTag)
    {
        const std::size_t num_rows = pGetLinearOperator(LhsMatrixTag)->Size1();
        const std::size_t num_cols = pGetLinearOperator(LhsMatrixTag)->Size2();
        const std::size_t size_dx = pGetVector(DxVectorTag)->size();
        const std::size_t size_rhs = pGetVector(RhsVectorTag)->size();
        return ((num_rows == num_cols) && (num_rows == size_dx) && (num_rows == size_rhs));
    }

    /**
     * @brief This method checks if the dimensions of the given linear system are consistent
     * @param LhsMatrixTag The tag of the system matrix to be checked (@see SparseMatrixTag)
     * @param RhsMatrixTag The tag of the residual matrix to be checked (@see DenseMatrixTag)
     * @param DxMatrixTag The tag of the solution matrix to be checked (@see DenseMatrixTag)
     * @return True if consistent, false otherwise
     */
    bool IsConsistent(
        SparseMatrixTag LhsMatrixTag,
        DenseMatrixTag RhsMatrixTag,
        DenseMatrixTag DxMatrixTag)
    {
        const std::size_t num_rows = pGetLinearOperator(LhsMatrixTag)->Size1();
        const std::size_t num_cols = pGetLinearOperator(LhsMatrixTag)->Size2();
        const std::size_t size_1_dx = pGetMatrix(DxMatrixTag)->size1();
        const std::size_t size_2_dx = pGetMatrix(DxMatrixTag)->size2();
        const std::size_t size_1_rhs = pGetMatrix(RhsMatrixTag)->size1();
        const std::size_t size_2_rhs = pGetMatrix(RhsMatrixTag)->size2();
        return ((num_rows == num_cols) && (num_rows == size_1_dx) && (num_rows == size_1_rhs) && (size_2_dx == size_2_rhs));
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
        return mpModelPart != nullptr && mpDofSet != nullptr;
    }

    ///@}
private:

    ///@name Member Variables
    ///@{

    const std::string mSystemName; // Name of the linear system

    ModelPart* mpModelPart = nullptr; // Model part of the linear system

    typename ModelPart::DofsArrayType* mpDofSet = nullptr; // Dofs of the linear system

    std::array<typename VectorType::Pointer, static_cast<std::size_t>(DenseVectorTag::NumberOfTags)> mVectors; // Vectors of the linear system

    std::array<DenseMatrixPointerType, static_cast<std::size_t>(DenseMatrixTag::NumberOfTags)> mDenseMatrices; // Dense matrices of the linear system

    std::array<typename LinearOperatorType::UniquePointer, static_cast<std::size_t>(SparseMatrixTag::NumberOfTags)> mLinearOperators; // Linear operators of the linear system

    ///@}
}; // Class LinearSystem

///@}
///@name Input and output
///@{

///@}
///@} addtogroup block

} // namespace Kratos.