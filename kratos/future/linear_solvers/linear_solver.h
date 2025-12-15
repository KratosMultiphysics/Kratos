//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//                   Ruben Zorrilla
//

#pragma once

// System includes

// External includes

// Project includes
#include "containers/csr_matrix.h"
#include "future/containers/linear_operator.h"
#include "containers/system_vector.h"
#include "includes/model_part.h"

namespace Kratos::Future
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class LinearSolver
 * @ingroup KratosCore
 * @brief Base class for all the linear solvers in Kratos.
 * @details This class define the general interface for the linear solvers in Kratos.
 * @tparam TSparseSpaceType which specify type of the unknowns, coefficients, sparse matrix, vector of unknowns, right hand side vector and their respective operators.
 * @tparam TDenseMatrixType which specify type of the matrices used as temporary matrices or multi solve unknowns and right hand sides and their operators.
 * @see SparseSpace
 * @see DenseSpace
 * @author Pooyan Dadvand
 * @author Riccardo Rossi
 * @author Ruben Zorrilla
 */
template<class TVectorType= SystemVector<>, typename... Ts>
class LinearSolver
{

static_assert(sizeof...(Ts) < 2, "Ts must be either a single type representing the sparse matrix type or empty to indicate a matrix-free operator.");

public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of LinearSolver
    KRATOS_CLASS_POINTER_DEFINITION(LinearSolver);

    /// Linear operator type definition
    using LinearOperatorType = LinearOperator<TVectorType, Ts...>;

    /// Linear operator pointer type definition
    using LinearOperatorPointerType = typename LinearOperatorType::Pointer;

    /// Sparse matrix type definition extracted from LinearOperator
    using MatrixPointerType = typename LinearOperatorType::SparseMatrixPointerType;

    /// Type definition for data
    using DataType = typename TVectorType::DataType;

    /// Type definition for index
    using IndexType = typename TVectorType::IndexType;

    /// Local system matrix type definition
    using DenseMatrixType = DenseMatrix<DataType>;

    /// Local system vector type definition
    using DenseTVectorType = DenseVector<DataType>;

    /// Boolean indicating whether a CSR matrix type is provided
    static constexpr bool kIsMatrixFree = LinearOperatorType::kIsMatrixFree;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    LinearSolver() {}

    /// Copy constructor.
    LinearSolver(const LinearSolver& Other) {}

    /// Destructor.
    virtual ~LinearSolver() = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    LinearSolver& operator=(const LinearSolver& Other)
    {
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function is designed to be called as few times as possible.
     * @details It creates the data structures that only depend on the connectivity of the matrix (and not on its coefficients) so that the memory can be allocated once and expensive operations can be done only when strictly  needed
     * @param pLinearOperator System matrix linear operator pointer.
     * @param rX Solution vector. it's also the initial guess for iterative linear solvers.
     * @param rB Right hand side vector.
     */
    virtual void Initialize(
        LinearOperatorPointerType pLinearOperator,
        TVectorType& rX,
        TVectorType& rB)
    {
    }

    /**
     * @brief This function is designed to be called every time the coefficients change in the system that is, normally at the beginning of each solve.
     * @details For example if we are implementing a direct solver, this is the place to do the factorization so that then the backward substitution can be performed effectively more than once
     * @param pLinearOperator System matrix linear operator pointer.
     * @param rX Solution vector. it's also the initial guess for iterative linear solvers.
     * @param rB Right hand side vector.
     */
    virtual void InitializeSolutionStep(
        LinearOperatorPointerType pLinearOperator,
        TVectorType& rX,
        TVectorType& rB)
    {
    }

    /**
     * @brief This function actually performs the solution work, eventually taking advantage of what was done before in the Initialize and InitializeSolutionStep functions.
     * @param pLinearOperator System matrix linear operator pointer.
     * @param rX Solution vector. it's also the initial guess for iterative linear solvers.
     * @param rB Right hand side vector.
     * @return @p true if the provided system was solved successfully satisfying the given requirements, @p false otherwise.
     */
    virtual bool PerformSolutionStep(
        LinearOperatorPointerType pLinearOperator,
        TVectorType& rX,
        TVectorType& rB)
    {
        KRATOS_ERROR << "Calling linear solver base class" << std::endl;
        return false;
    }

    /**
     * @brief This function is designed to be called at the end of the solve step.
     * @details for example this is the place to remove any data that we do not want to save for later
     * @param pLinearOperator System matrix linear operator pointer.
     * @param rX Solution vector. it's also the initial guess for iterative linear solvers.
     * @param rB Right hand side vector.
     */
    virtual void FinalizeSolutionStep(
        LinearOperatorPointerType pLinearOperator,
        TVectorType& rX,
        TVectorType& rB)
    {
    }

    /**
     * @brief This function is designed to clean up all internal data in the solver.
     * @details Clear is designed to leave the solver object as if newly created. After a clear a new Initialize is needed
     */
    virtual void Clear()
    {
    }

    /**
     * @brief Solve method with linear operator as input.
     * @details Solves the linear system Ax=b and puts the result on SystemVector& rX. rX is also th initial guess for iterative methods.
     * @param rLinearOperator System matrix linear operator
     * @param rX Solution vector
     * @param rB Right hand side vector
     * @return true if the system was solved successfully, false otherwise
     */
    virtual bool Solve(
        LinearOperatorPointerType pLinearOperator,
        TVectorType& rX,
        TVectorType& rB)
    {
        this->InitializeSolutionStep(pLinearOperator, rX, rB);
        const auto status = this->PerformSolutionStep(pLinearOperator, rX, rB);
        this->FinalizeSolutionStep(pLinearOperator, rX, rB);
        this->Clear();
        return status;
    }

    /**
     * @brief Solve method overload taking a sparse matrix as input.
     * @details Solves the linear system Ax=b and puts the result on SystemVector& rX. rX is also the initial guess for iterative methods.
     * Note that this method internally creates a LinearOperator from the provided sparse matrix and calls the corresponding Solve method.
     * @param pA Matrix type pointer representing the system
     * @param rX Solution vector
     * @param rB Right hand side vector
     * @return true if the system was solved successfully, false otherwise
     */
    bool Solve(
        MatrixPointerType pA,
        TVectorType& rX,
        TVectorType& rB)
    {
        auto p_linear_operator = Kratos::make_shared<LinearOperatorType>(*pA);
        return Solve(p_linear_operator, rX, rB);
    }

    /**
     * @brief Multi solve method for solving a set of linear systems with same coefficient matrix.
     * @details Solves the linear system Ax=b and puts the result on SystemVector& rX. rX is also th initial guess for iterative methods.
     * @param rLinearOperator System matrix linear operator
     * @param rX Solution vector
     * @param rB Right hand side vector
     * @return true if the system was solved successfully, false otherwise
     */
    virtual bool Solve(
        LinearOperatorPointerType rLinearOperator,
        DenseMatrixType& rX,
        DenseMatrixType& rB)
    {
        KRATOS_ERROR << "Calling linear solver base class" << std::endl;
        return false;
    }

    /**
     * @brief Multi solve method for solving a set of linear systems with same coefficient matrix.
     * @details Solves the linear system Ax=b and puts the result on SystemVector& rX. rX is also the initial guess for iterative methods.
     * Note that this method internally creates a LinearOperator from the provided sparse matrix and calls the corresponding Solve method.
     * @param pA Pointer to the sparse matrix representing the system
     * @param rX Solution vector
     * @param rB Right hand side vector
     * @return true if the system was solved successfully, false otherwise
     */
    bool Solve(
        MatrixPointerType pA,
        DenseMatrixType& rX,
        DenseMatrixType& rB)
    {
        auto p_linear_operator = Kratos::make_shared<LinearOperatorType>(*pA);
        return Solve(p_linear_operator, rX, rB);
    }

    /**
     * @brief Eigenvalue and eigenvector solve method for derived eigensolvers
     * @param K The stiffness matrix linear operator
     * @param M The mass matrix linear operator
     * @param Eigenvalues The vector containing the eigen values
     * @param Eigenvectors The matrix containing the eigen vectors
     */
    virtual void Solve(
        LinearOperatorPointerType K,
        LinearOperatorPointerType M,
        DenseTVectorType& Eigenvalues,
        DenseMatrixType& Eigenvectors)
    {
        KRATOS_ERROR << "Calling linear solver base class" << std::endl;
    }

    /**
     * @brief Eigenvalue and eigenvector solve method for derived eigensolvers
     * Note that this method internally creates a LinearOperator from the provided stiffness and mass matrices and calls the corresponding Solve method.
     * @param pK Pointer to the stiffness matrix
     * @param pM Pointer to the mass matrix
     * @param Eigenvalues The vector containing the eigen values
     * @param Eigenvectors The matrix containing the eigen vectors
     */
    void Solve(
        MatrixPointerType pK,
        MatrixPointerType pM,
        DenseTVectorType& Eigenvalues,
        DenseMatrixType& Eigenvectors)
    {
        auto p_linear_operator_K = Kratos::make_shared<LinearOperatorType>(*pK);
        auto p_linear_operator_M = Kratos::make_shared<LinearOperatorType>(*pM);
        Solve(p_linear_operator_K, p_linear_operator_M, Eigenvalues, Eigenvectors);
    }

    /**
     * @brief Checks if additional physical data is needed by the solver.
     * @details Some solvers may require a minimum degree of knowledge of the structure of the matrix.
     * For instance, when solving a mixed u-p problem, it is important to identify the row associated with v and p.
     * Another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers,
     * which require knowledge of the spatial position of the nodes associated with a given degree of freedom (DOF).
     * @return True if additional physical data is needed, false otherwise.
     */
    virtual bool AdditionalPhysicalDataIsNeeded()
    {
        return false;
    }

    /**
     * @brief Provides additional physical data required by the solver.
     * @details Some solvers may require a minimum degree of knowledge of the structure of the matrix.
     * For example, when solving a mixed u-p problem, it is important to identify the row associated with v and p.
     * Another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers,
     * which require knowledge of the spatial position of the nodes associated with a given degree of freedom (DOF).
     * This function provides the opportunity to provide such data if needed.
     * @param pLinearOperator Pointer to the sparse matrix linear operator.
     * @param rX The solution vector.
     * @param rB The right-hand side vector.
     * @param rDoFSet The set of degrees of freedom.
     * @param rModelPart The model part.
     */
    virtual void ProvideAdditionalData(
        LinearOperatorPointerType pLinearOperator,
        TVectorType& rX,
        TVectorType& rB,
        typename ModelPart::DofsArrayType& rDoFSet,
        ModelPart& rModelPart)
    {
    }

    /**
     * @brief Provides additional physical data required by the solver.
     * @details Some solvers may require a minimum degree of knowledge of the structure of the matrix.
     * For example, when solving a mixed u-p problem, it is important to identify the row associated with v and p.
     * Another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers,
     * which require knowledge of the spatial position of the nodes associated with a given degree of freedom (DOF).
     * This function provides the opportunity to provide such data if needed.
     * Note that this method internally creates a LinearOperator from the provided sparse matrix and calls the corresponding Solve method.
     * @param pA Pointer to the sparse matrix.
     * @param rX The solution vector.
     * @param rB The right-hand side vector.
     * @param rDoFSet The set of degrees of freedom.
     * @param rModelPart The model part.
     */
    void ProvideAdditionalData(
        MatrixPointerType pA,
        TVectorType& rX,
        TVectorType& rB,
        typename ModelPart::DofsArrayType& rDoFSet,
        ModelPart& rModelPart)
    {
        auto p_linear_operator = Kratos::make_shared<LinearOperatorType>(*pA);
        return ProvideAdditionalData(p_linear_operator, rX, rB, rDoFSet, rModelPart);
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief This method allows to set the tolerance in the linear solver
     * @param NewTolerance The new tolerance set
     */
    virtual void SetTolerance(double NewTolerance)
    {
        KRATOS_WARNING("LinearSolver") << "Accessed base function \"SetTolerance\". This does nothing !" << std::endl;
    }

    /**
     * @brief This method allows to get the tolerance in the linear solver
     * @return The tolerance
     */
    virtual double GetTolerance()
    {
        KRATOS_WARNING("LinearSolver") << "Accessed base function \"GetTolerance\". No tolerance defined, returning 0 !" << std::endl ;
        return 0;
    }

    /**
     * @brief Virtual function to get the number of iterations.
     * @details This function returns the number of iterations performed by the linear solver. As this is a base function, it returns 0 by default and issues a warning.
     * @return The number of iterations performed.
     */
    virtual IndexType GetIterationsNumber()
    {
        KRATOS_WARNING("LinearSolver") << "Accessed base function \"GetIterationsNumber\", returning 0 !" << std::endl;

        return 0;
    }

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief This method checks if the dimensions of the system of equations are consistent
     * @param pLinearOperator Pointer to the sparse matrix linear operator
     * @param rX The vector containing the unknowns
     * @param rB The RHS of the system of equations
     * @return True if consistent, false otherwise
     */
    virtual bool IsConsistent(
        LinearOperatorPointerType pLinearOperator,
        TVectorType& rX,
        TVectorType& rB)
    {
        const std::size_t num_rows = pLinearOperator->NumRows();
        const std::size_t num_cols = pLinearOperator->NumCols();
        const std::size_t size_x = rX.size();
        const std::size_t size_b = rB.size();

        return ((num_rows ==  num_cols) && (num_rows ==  size_x) && (num_rows ==  size_b));
    }

    /**
     * @brief This method checks if the dimensions of the system of equations are consistent (dense matrix for RHS and unknowns version)
     * @param pLinearOperator Pointer to the sparse matrix linear operator
     * @param rX The matrix containing the unknowns
     * @param rB The matrix containing the RHSs of the system of equations
     * @return True if consistent, false otherwise
     */
    virtual bool IsConsistent(
        LinearOperatorPointerType pLinearOperator,
        DenseMatrixType& rX,
        DenseMatrixType& rB)
    {
        const std::size_t num_rows = pLinearOperator->NumRows();
        const std::size_t num_cols = pLinearOperator->NumCols();
        const std::size_t size_1_x = rX.size1();
        const std::size_t size_1_b = rB.size1();
        const std::size_t size_2_x = rX.size2();
        const std::size_t size_2_b = rB.size2();

        return ((num_rows ==  num_cols) && (num_rows ==  size_1_x) && (num_rows ==  size_1_b) && (size_2_x == size_2_b));
    }

    /**
     * @brief This method checks if the dimensions of the system of equations are not consistent
     * @param pLinearOperator Pointer to the sparse matrix linear operator
     * @param rX The vector containing the unknowns
     * @param rB The RHS of the system of equations
     * @return False if consistent, true otherwise
     */
    virtual bool IsNotConsistent(
        LinearOperatorPointerType pLinearOperator,
        TVectorType& rX,
        TVectorType& rB)
    {
        return !IsConsistent(pLinearOperator, rX, rB);
    }

    /**
     * @brief This method checks if the dimensions of the system of equations are not consistent
     * @param pLinearOperator Pointer to the sparse matrix linear operator
     * @param rX The matrix containing the unknowns
     * @param rB The matrix containing the RHSs of the system of equations
     * @return False if consistent, true otherwise
     */
    virtual bool IsNotConsistent(
        LinearOperatorPointerType pLinearOperator,
        DenseMatrixType& rX,
        DenseMatrixType& rB)
    {
        return !IsConsistent(pLinearOperator, rX, rB);
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "New linear solver";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "New linear solver";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}
};

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template<class TVectorType, typename... Ts>
inline std::istream& operator >> (
    std::istream& IStream,
    LinearSolver<TVectorType, Ts...>& rThis)
{
    return IStream;
}

/// output stream function
template<class TVectorType, typename... Ts>
inline std::ostream& operator << (
    std::ostream& rOStream,
    const LinearSolver<TVectorType, Ts...>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos::Future.
