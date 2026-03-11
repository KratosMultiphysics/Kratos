//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#pragma once

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "future/linear_solvers/iterative_solver.h"

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
 * @class CGSolver
 * @ingroup KratosCore
 * @brief Preconditioned Conjugate Gradient (PCG) linear solver.
 * @details This class implements the Preconditioned Conjugate Gradient (PCG) iterative method for solving linear systems.
 *          It is derived from the IterativeSolver base class.
 * @tparam TLinearAlgebra Type of the linear algebra used (e.g. CSR space, dense space).
 */
template<class TLinearAlgebra>
class CGSolver : public Future::IterativeSolver<TLinearAlgebra>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CGSolver
    KRATOS_CLASS_POINTER_DEFINITION(CGSolver);

    /// Base type iterative solver definition
    using BaseType = Future::IterativeSolver<TLinearAlgebra>;

    /// Sparse matrix type definition from linear algebra traits
    using MatrixType = typename TLinearAlgebra::MatrixType;

    /// Dense vector type definition from linear algebra traits
    using VectorType = typename TLinearAlgebra::VectorType;

    /// Dense matrix type definition from linear algebra traits
    using DenseMatrixType = typename TLinearAlgebra::DenseMatrixType;

    /// Linear system container type definition
    using LinearSystemType = LinearSystem<TLinearAlgebra>;

    /// Linear operator type definition
    using LinearOperatorType = LinearOperator<TLinearAlgebra>;

    /// Preconditions pointer type definition
    using PreconditionerPointerType = typename Preconditioner<TLinearAlgebra>::Pointer;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CGSolver(Parameters Settings = Parameters(R"({})"))
        : BaseType(Settings)
    {
        std::cout << "CG constructor" << std::endl;
        KRATOS_WATCH(Settings)
        // Validate and assign default parameters
        Settings.ValidateAndAssignDefaults(BaseType::GetDefaultParameters());
        KRATOS_WATCH(Settings)


    }

    /// Constructor with preconditioner.
    CGSolver(
        Parameters Settings,
        PreconditionerPointerType pPreconditioner)
        : BaseType(Settings, pPreconditioner)
    {
    }

    /// Copy constructor.
    CGSolver(const CGSolver& Other) = delete;

    /// Destructor.
    ~CGSolver() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    CGSolver& operator=(const CGSolver& Other) = delete;

    ///@}
    ///@name Operations
    ///@{

    bool PerformSolutionStep(LinearSystemType& rLinearSystem) override
    {
        bool is_solved;
        if (!this->mMultipleSolve) {
            // Get sparse matrix and dense vector tags from strings
            const auto dx_tag = LinearSystemTags::DenseVectorTagFromString(BaseType::mDxTagString);
            const auto rhs_tag = LinearSystemTags::DenseVectorTagFromString(BaseType::mRhsTagString);
            const auto lhs_tag = LinearSystemTags::SparseMatrixTagFromString(BaseType::mLhsTagString);

            // Check if the linear system is consistent
            if(!rLinearSystem.IsConsistent(lhs_tag, rhs_tag, dx_tag)) {
                KRATOS_WARNING("CGSolver") << "Linear system is not consistent. PerformSolutionStep cannot be performed." << std::endl;
                return false;
            }

            // Get arrays from linear system container
            auto& r_dx = *(rLinearSystem.pGetVector(dx_tag));
            auto& r_rhs = *(rLinearSystem.pGetVector(rhs_tag));
            const auto& rp_lhs_lin_op = rLinearSystem.pGetLinearOperator(lhs_tag);

            // Call CG implementation
            is_solved = IterativeSolve(rp_lhs_lin_op, r_rhs, r_dx);
            KRATOS_WARNING_IF("CGSolver", !is_solved) << "Non converged linear solution. ["<< BaseType::GetResidualNorm() / BaseType::mBNorm << " > "<<  BaseType::GetTolerance() << "]" << std::endl;

        } else {
            // Get sparse matrix and dense matrices tags from strings
            const auto dx_tag = LinearSystemTags::DenseMatrixTagFromString(BaseType::mDxTagString);
            const auto rhs_tag = LinearSystemTags::DenseMatrixTagFromString(BaseType::mRhsTagString);
            const auto lhs_tag = LinearSystemTags::SparseMatrixTagFromString(BaseType::mLhsTagString);

            // Get arrays from linear system container
            auto& r_dx = *(rLinearSystem.pGetMatrix(dx_tag));
            const auto& r_rhs = *(rLinearSystem.pGetMatrix(rhs_tag));
            const auto& rp_lhs_lin_op = rLinearSystem.pGetLinearOperator(lhs_tag);

            // Check if the linear system is consistent
            if(!rLinearSystem.IsConsistent(lhs_tag, rhs_tag, dx_tag)) {
                KRATOS_WARNING("CGSolver") << "Linear system is not consistent. PerformSolutionStep cannot be performed." << std::endl;
                return false;
            }

            // Perform the columnwise solve
            is_solved = true;
            VectorType aux_dx(r_rhs.size1()); // Auxiliary update vector
            VectorType aux_rhs(r_rhs.size1()); // Auxiliary residual vector
            const std::size_t n_problems = r_rhs.size2();
            for (std::size_t i = 0; i < n_problems; ++i) {
                // Get current residual column
                IndexPartition<IndexType>(aux_rhs.size()).for_each([i, &r_rhs, &aux_rhs](const auto Index) {
                    aux_rhs[Index] = r_rhs(i, Index);
                });

                // Call CG implementation
                const bool aux_is_solved = IterativeSolve(rp_lhs_lin_op, aux_rhs, aux_dx);
                KRATOS_WARNING_IF("CGSolver", !aux_is_solved) << "Non converged linear solution for column " << i << ". ["<< BaseType::GetResidualNorm() / BaseType::mBNorm << " > "<<  BaseType::GetTolerance() << "]" << std::endl;
                is_solved &= aux_is_solved;

                // Set current solution update
                IndexPartition<IndexType>(aux_dx.size()).for_each([i, &r_dx, &aux_dx](const auto Index) {
                    r_dx(Index, i) = aux_dx[Index];
                });
            }
        }

        return is_solved;
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Conjugate gradient linear solver with " << BaseType::GetPreconditioner()->Info();
        return  buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
    }

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

    bool IterativeSolve(
        const typename LinearOperatorType::UniquePointer& rpLinearOperator,
        const VectorType& rB,
        VectorType& rX)
    {
        // Initialize variables
        const int size = rX.size(); // Size of the system
        BaseType::mBNorm = rB.Norm(); // Note that this will be used to compute the relative residual (@see BaseType::IterationNeeded)
        BaseType::mIterationsNumber = 0;

        // r = rB - A * rX
        VectorType r(size);
        r.SetValue(0.0);
        rpLinearOperator->SpMV(rX, r);
        r *= -1.0;
        r += rB;
        BaseType::mResidualNorm = r.Norm();

        // z = M^{-1} * r
        VectorType z(size);
        z.SetValue(0.0);
        this->GetPreconditioner()->Apply(r, z);

        // Initialize iteration variables
        VectorType p(z);
        VectorType q(size);
        q.SetValue(0.0);

        double rho_0 = r.Dot(z);
        double rho_1 = rho_0;

        // Early return condition for zero residual
        if (std::abs(rho_0) < std::numeric_limits<double>::epsilon()) {
            return true;
        }

        // Main loop
        while(BaseType::IterationNeeded() && (std::abs(rho_0) > std::numeric_limits<double>::epsilon())) {

            // q = A * p
            q.SetValue(0.0);
            rpLinearOperator->SpMV(p, q);

            // pq = p^{T} * A * p
            const double pq = p.Dot(q);
            if (std::abs(pq) < std::numeric_limits<double>::epsilon()) {
                break;
            }

            // alpha: step length along the conjugate direction p
            const double alpha = rho_0 / pq;

            // dx = dx + alpha * p
            rX.Add(alpha, p);

            // r = r - alpha * q
            r.Add(-alpha, q);
            BaseType::mResidualNorm = r.Norm(); // Update residual norm for convergence check (@see BaseType::IterationNeeded)

            // z = M^{-1} * r
            z.SetValue(0.0);
            this->GetPreconditioner()->Apply(r, z);

            // rho_1 = r^{T} * z
            rho_1 = r.Dot(z);

            // beta: coefficient to build the next conjugate search direction
            const double beta = (rho_1 / rho_0);

            // p = z + beta * p
            p *= beta;
            p.Add(1.0, z);

            // Next iteration update
            rho_0 = rho_1;
            BaseType::mIterationsNumber++; // Update iterations counter for convergence check (@see BaseType::IterationNeeded)
        }

        return BaseType::IsConverged();
    }

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
}; // Class CGSolver

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
    CGSolver<TLinearAlgebra>& rThis)
{
    return IStream;
}

/// output stream function
template<class TLinearAlgebra>
inline std::ostream& operator << (
    std::ostream& OStream,
    const CGSolver<TLinearAlgebra>& rThis)
{
    rThis.PrintInfo(OStream);
    OStream << std::endl;
    rThis.PrintData(OStream);

    return OStream;
}

///@}

}  // namespace Kratos.



