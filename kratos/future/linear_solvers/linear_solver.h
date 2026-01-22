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
#include "containers/csr_matrix.h"
#include "future/containers/define_linear_algebra_serial.h"
#include "future/containers/linear_system.h"
#include "future/containers/eigenvalue_system.h"
#include "includes/model_part.h"
#include "future/linear_operators/sparse_matrix_linear_operator.h"

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
template<class TLinearAlgebra>
class LinearSolver
{

public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of LinearSolver
    KRATOS_CLASS_POINTER_DEFINITION(LinearSolver);

    /// Linear system type definition
    using LinearSystemType = LinearSystem<TLinearAlgebra>;

    /// Eigenvalue system type definition
    using EigenvalueSystemType = EigenvalueSystem<TLinearAlgebra>;

    /// Linear operator pointer type definition
    using LinearOperatorPointerType = typename LinearOperator<TLinearAlgebra>::Pointer;

    /// Vector type definition from linear algebra template parameter
    using VectorType = typename TLinearAlgebra::VectorType;

    /// Dense matrix type definition from linear algebra template parameter used for multi-solve and eigenvector storage
    using DenseMatrixType = typename TLinearAlgebra::DenseMatrixType;

    /// Type definition for data
    using DataType = typename VectorType::DataType;

    /// Type definition for index
    using IndexType = typename VectorType::IndexType;

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
     * @details It creates the data structures that only depend on the connectivity of the matrix (and not on its coefficients) so that the memory can be allocated once and expensive operations can be done only when strictly needed
     * @param rLinearSystem The linear system to be solved
     */
    virtual void Initialize(LinearSystemType& rLinearSystem)
    {
    }

    /**
     * @brief This function is designed to be called as few times as possible.
     * @details It creates the data structures that only depend on the connectivity of the matrix (and not on its coefficients) so that the memory can be allocated once and expensive operations can be done only when strictly needed
     * @param rEigenvalueSystem The eigenvlaue system to be solved
     */
    virtual void Initialize(EigenvalueSystemType& rEigenvalueSystem)
    {
    }

    //TODO:
    // virtual void Initialize(MultiRHSLinearSystemType>& rLinearSystem)
    // {
    // }

    /**
     * @brief This function is designed to be called every time the coefficients change in the system that is, normally at the beginning of each solve.
     * @details For example if we are implementing a direct solver, this is the place to do the factorization so that then the backward substitution can be performed effectively more than once
     * @param rLinearSystem The linear system to be solved
     */
    virtual void InitializeSolutionStep(LinearSystemType& rLinearSystem)
    {
    }
    /**
     * @brief This function is designed to be called every time the coefficients change in the system that is, normally at the beginning of each solve.
     * @details For example if we are implementing a direct solver, this is the place to do the factorization so that then the backward substitution can be performed effectively more than once
     * @param rEigenvalueSystem The eigenvalue system to be solved
     */
    virtual void InitializeSolutionStep(EigenvalueSystemType& rEigenvalueSystem)
    {
    }

    /**
     * @brief This function actually performs the solution work, eventually taking advantage of what was done before in the Initialize and InitializeSolutionStep functions.
     * @param rLinearSystem The linear system to be solved.
     * @return @p true if the provided system was solved successfully satisfying the given requirements, @p false otherwise.
     */
    virtual bool PerformSolutionStep(LinearSystemType& rLinearSystem)
    {
        KRATOS_ERROR << "Calling linear solver base class." << std::endl;
        return false;
    }

    /**
     * @brief This function actually performs the solution work, eventually taking advantage of what was done before in the Initialize and InitializeSolutionStep functions.
     * @param rLinearSystem The linear system to be solved.
     * @return @p true if the provided system was solved successfully satisfying the given requirements, @p false otherwise.
     */
    virtual bool PerformSolutionStep(EigenvalueSystemType& rEigenvalueSystem)
    {
        KRATOS_ERROR << "Calling linear solver base class." << std::endl;
        return false;
    }

    /**
     * @brief This function is designed to be called at the end of the solve step.
     * @details for example this is the place to remove any data that we do not want to save for later
     * @param rLinearSystem The linear system to be solved.
     */
    virtual void FinalizeSolutionStep(LinearSystemType& rLinearSystem)
    {
    }

    /**
     * @brief This function is designed to be called at the end of the solve step.
     * @details for example this is the place to remove any data that we do not want to save for later
     * @param rEigenvalueSystem The linear system to be solved.
     */
    virtual void FinalizeSolutionStep(EigenvalueSystemType& rEigenvalueSystem)
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
        VectorType& rX,
        VectorType& rB,
        typename ModelPart::DofsArrayType& rDoFSet,
        ModelPart& rModelPart)
    {
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

    //TODO: move this to the multivector linear system once we implement it!!!!
    // /**
    //  * @brief This method checks if the dimensions of the system of equations are consistent (dense matrix for RHS and unknowns version)
    //  * @param pLinearOperator Pointer to the sparse matrix linear operator
    //  * @param rX The matrix containing the unknowns
    //  * @param rB The matrix containing the RHSs of the system of equations
    //  * @return True if consistent, false otherwise
    //  */
    // virtual bool IsConsistent(
    //     LinearOperatorPointerType pLinearOperator,
    //     DenseMatrixType& rX,
    //     DenseMatrixType& rB)
    // {
    //     const std::size_t num_rows = pLinearOperator->NumRows();
    //     const std::size_t num_cols = pLinearOperator->NumCols();
    //     const std::size_t size_1_x = rX.size1();
    //     const std::size_t size_1_b = rB.size1();
    //     const std::size_t size_2_x = rX.size2();
    //     const std::size_t size_2_b = rB.size2();

    //     return ((num_rows ==  num_cols) && (num_rows ==  size_1_x) && (num_rows ==  size_1_b) && (size_2_x == size_2_b));
    // }

    // /**
    //  * @brief This method checks if the dimensions of the system of equations are not consistent
    //  * @param pLinearOperator Pointer to the sparse matrix linear operator
    //  * @param rX The vector containing the unknowns
    //  * @param rB The RHS of the system of equations
    //  * @return False if consistent, true otherwise
    //  */
    // virtual bool IsNotConsistent(
    //     LinearOperatorPointerType pLinearOperator,
    //     VectorType& rX,
    //     VectorType& rB)
    // {
    //     return !IsConsistent(pLinearOperator, rX, rB);
    // }

    // /**
    //  * @brief This method checks if the dimensions of the system of equations are not consistent
    //  * @param pLinearOperator Pointer to the sparse matrix linear operator
    //  * @param rX The matrix containing the unknowns
    //  * @param rB The matrix containing the RHSs of the system of equations
    //  * @return False if consistent, true otherwise
    //  */
    // virtual bool IsNotConsistent(
    //     LinearOperatorPointerType pLinearOperator,
    //     DenseMatrixType& rX,
    //     DenseMatrixType& rB)
    // {
    //     return !IsConsistent(pLinearOperator, rX, rB);
    // }

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
template<class TLinearAlgebra>
inline std::istream& operator >> (
    std::istream& IStream,
    LinearSolver<TLinearAlgebra>& rThis)
{
    return IStream;
}

/// output stream function
template<class TLinearAlgebra>
inline std::ostream& operator << (
    std::ostream& rOStream,
    const LinearSolver<TLinearAlgebra>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos::Future.
