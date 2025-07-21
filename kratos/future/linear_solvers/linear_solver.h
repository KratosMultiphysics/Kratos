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
template<class TMatrixType = CsrMatrix<>, class TVectorType= SystemVector<>>
class LinearSolver
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of LinearSolver
    KRATOS_CLASS_POINTER_DEFINITION(LinearSolver);

    /// Type definition for sparse matrix
    using SparseMatrixType = TMatrixType;

    /// Type definition for pointer to sparse matrix
    using SparseMatrixPointerType = typename SparseMatrixType::Pointer;

    /// Type definition for vector
    using VectorType = TVectorType;

    /// Type definition for pointer to vector
    using VectorPointerType = typename VectorType::Pointer;

    /// Type definition for data
    using DataType = typename SparseMatrixType::DataType;

    /// Type definition for size
    using SizeType = typename SparseMatrixType::DataType;

    /// Type definition for index
    using IndexType = typename SparseMatrixType::IndexType;

    /// Local system matrix type definition
    using DenseMatrixType = DenseMatrix<DataType>;

    /// Local system vector type definition
    using DenseVectorType = DenseVector<DataType>;

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
     * @param rA. System matrix
     * @param rX. Solution vector. it's also the initial guess for iterative linear solvers.
     * @param rB. Right hand side vector.
     */
    virtual void Initialize(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
    }

    /**
     * @brief This function is designed to be called every time the coefficients change in the system that is, normally at the beginning of each solve.
     * @details For example if we are implementing a direct solver, this is the place to do the factorization so that then the backward substitution can be performed effectively more than once
     * @param rA. System matrix
     * @param rX. Solution vector. it's also the initial guess for iterative linear solvers.
     * @param rB. Right hand side vector.
     */
    virtual void InitializeSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
    }

    /**
     * @brief This function actually performs the solution work, eventually taking advantage of what was done before in the Initialize and InitializeSolutionStep functions.
     * @param rA. System matrix
     * @param rX. Solution vector. it's also the initial guess for iterative linear solvers.
     * @param rB. Right hand side vector.
     * @return @p true if the provided system was solved successfully satisfying the given requirements, @p false otherwise.
     */
    virtual bool PerformSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
        KRATOS_ERROR << "Calling linear solver base class" << std::endl;
        return false;
    }

    /**
     * @brief This function is designed to be called at the end of the solve step.
     * @details for example this is the place to remove any data that we do not want to save for later
     * @param rA. System matrix
     * @param rX. Solution vector. it's also the initial guess for iterative linear solvers.
     * @param rB. Right hand side vector.
     */
    virtual void FinalizeSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
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
     * @brief Normal solve method.
     * @details Solves the linear system Ax=b and puts the result on SystemVector& rX. rX is also th initial guess for iterative methods.
     * @param rA. System matrix
     * @param rX. Solution vector. it's also the initial
    guess for iterative linear solvers.
     * @param rB. Right hand side vector.
     */
    virtual bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
        KRATOS_ERROR << "Calling linear solver base class" << std::endl;
        return false;
    }

    /**
     * @brief Multi solve method for solving a set of linear systems with same coefficient matrix.
     * @details Solves the linear system Ax=b and puts the result on SystemVector& rX. rX is also th initial guess for iterative methods.
     * @param rA. System matrix
     * @param rX. Solution vector. it's also the initial
    guess for iterative linear solvers.
     * @param rB. Right hand side vector.
     */
    virtual bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
    {
        KRATOS_ERROR << "Calling linear solver base class" << std::endl;
        return false;
    }

    /**
     * @brief Eigenvalue and eigenvector solve method for derived eigensolvers
     * @param K The stiffness matrix
     * @param M The mass matrix
     * @param Eigenvalues The vector containing the eigen values
     * @param Eigenvectors The matrix containing the eigen vectors
     */
    virtual  void Solve(SparseMatrixType& K,
                        SparseMatrixType& M,
                        DenseVectorType& Eigenvalues,
                        DenseMatrixType& Eigenvectors)
    {
        KRATOS_ERROR << "Calling linear solver base class" << std::endl;
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
     * @param rA The sparse matrix.
     * @param rX The solution vector.
     * @param rB The right-hand side vector.
     * @param rDoFSet The set of degrees of freedom.
     * @param rModelPart The model part.
     */
    virtual void ProvideAdditionalData(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB,
        typename ModelPart::DofsArrayType& rDoFSet,
        ModelPart& rModelPart
        )
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

    /**
     * @brief This method checks if the dimensions of the system of equations are consistent
     * @param rA The LHS of the system of equations
     * @param rX The vector containing the unknowns
     * @param rB The RHS of the system of equations
     * @return True if consistent, false otherwise
     */
    virtual bool IsConsistent(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB)
    {
        const SizeType size = rA.size1();
        const SizeType size_a = rA.size2();
        const SizeType size_x = rX.size();
        const SizeType size_b = rB.size();

        return ((size ==  size_a) &&
                (size ==  size_x) &&
                (size ==  size_b));
    }

    /**
     * @brief This method checks if the dimensions of the system of equations are consistent (dense matrix for RHS and unknowns version)
     * @param rA The LHS of the system of equations
     * @param rX The matrix containing the unknowns
     * @param rB The matrix containing the RHSs of the system of equations
     * @return True if consistent, false otherwise
     */
    virtual bool IsConsistent(
        SparseMatrixType& rA,
        DenseMatrixType& rX,
        DenseMatrixType& rB)
    {
        const SizeType size = rA.size1();
        const SizeType size_a = rA.size2();
        const SizeType size_1_x = rX.size1();
        const SizeType size_1_b = rB.size1();
        const SizeType size_2_x = rX.size2();
        const SizeType size_2_b = rB.size2();

        return ((size ==  size_a) &&
                (size ==  size_1_x) &&
                (size ==  size_1_b) &&
                (size_2_x == size_2_b));
    }

    /**
     * @brief This method checks if the dimensions of the system of equations are not consistent
     * @param rA The LHS of the system of equations
     * @param rX The vector containing the unknowns
     * @param rB The RHS of the system of equations
     * @return False if consistent, true otherwise
     */
    virtual bool IsNotConsistent(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB)
    {
        return (!IsConsistent(rA, rX, rB));
    }

    /**
     * @brief This method checks if the dimensions of the system of equations are not consistent
     * @param rA The LHS of the system of equations
     * @param rX The matrix containing the unknowns
     * @param rB The matrix containing the RHSs of the system of equations
     * @return False if consistent, true otherwise
     */
    virtual bool IsNotConsistent(
        SparseMatrixType& rA,
        DenseMatrixType& rX,
        DenseMatrixType& rB)
    {
        return (!IsConsistent(rA, rX, rB));
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
template<class TMatrixType, class TVectorType>
inline std::istream& operator >> (
    std::istream& IStream,
    LinearSolver<TMatrixType, TVectorType>& rThis)
{
    return IStream;
}

/// output stream function
template<class TMatrixType, class TVectorType>
inline std::ostream& operator << (
    std::ostream& rOStream,
    const LinearSolver<TMatrixType, TMatrixType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos::Future.
