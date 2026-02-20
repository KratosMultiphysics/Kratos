//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
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

/// Short class definition.
/** Detail class definition.
*/
template<class TLinearAlgebra>
class CGSolver : public IterativeSolver<TLinearAlgebra>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CGSolver
    KRATOS_CLASS_POINTER_DEFINITION(CGSolver);

    /// Base type iterative solver definition
    using BaseType = IterativeSolver<TLinearAlgebra>;

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
        bool is_solved = false;
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

            // Apply preconditioner
            auto p_preconditioner = this->GetPreconditioner();
            p_preconditioner->ApplyInverseRight(r_dx);
            p_preconditioner->ApplyLeft(r_rhs);

            is_solved = IterativeSolve(rp_lhs_lin_op, r_rhs, r_dx);

            KRATOS_WARNING_IF("CGSolver", !is_solved) << "Non converged linear solution. ["<< BaseType::GetResidualNorm() / BaseType::mBNorm << " > "<<  BaseType::GetTolerance() << "]" << std::endl;

            // BaseType::GetPreconditioner()->Finalize(rX);


        } else {
            // Get sparse matrix and dense matrices tags from strings
            const auto dx_tag = LinearSystemTags::DenseMatrixTagFromString(BaseType::mDxTagString);
            const auto rhs_tag = LinearSystemTags::DenseMatrixTagFromString(BaseType::mRhsTagString);
            const auto lhs_tag = LinearSystemTags::SparseMatrixTagFromString(BaseType::mLhsTagString);

            // Check if the linear system is consistent
            if(!rLinearSystem.IsConsistent(lhs_tag, rhs_tag, dx_tag)) {
                KRATOS_WARNING("CGSolver") << "Linear system is not consistent. PerformSolutionStep cannot be performed." << std::endl;
                return false;
            }

            KRATOS_ERROR << "Multiple solve not implemented for CG solver" << std::endl;
        }

        return is_solved;
    }


//     /** Normal solve method.
//     Solves the linear system Ax=b and puts the result on SystemVector& rX.
//     rX is also th initial guess for iterative methods.
//     @param rA. System matrix
//     @param rX. Solution vector. it's also the initial
//     guess for iterative linear solvers.
//     @param rB. Right hand side vector.
//     */
//     bool Solve(MatrixType& rA, VectorType& rX, VectorType& rB) override
//     {
//         if(this->IsNotConsistent(rA, rX, rB))
//             return false;

// // 	  GetTimeTable()->Start(Info());

//         BaseType::GetPreconditioner()->Initialize(rA,rX,rB);
//         BaseType::GetPreconditioner()->ApplyInverseRight(rX);
//         BaseType::GetPreconditioner()->ApplyLeft(rB);

//         bool is_solved = IterativeSolve(rA,rX,rB);

//         KRATOS_WARNING_IF("CG Linear Solver", !is_solved)<<"Non converged linear solution. ["<< BaseType::GetResidualNorm()/BaseType::mBNorm << " > "<<  BaseType::GetTolerance() << "]" << std::endl;

//         BaseType::GetPreconditioner()->Finalize(rX);

// // 	  GetTimeTable()->Stop(Info());

//         return is_solved;
//     }


//     /** Multi solve method for solving a set of linear systems with same coefficient matrix.
//     Solves the linear system Ax=b and puts the result on SystemVector& rX.
//     rX is also th initial guess for iterative methods.
//     @param rA. System matrix
//     @param rX. Solution vector. it's also the initial
//     guess for iterative linear solvers.
//     @param rB. Right hand side vector.
//     */
//     bool Solve(MatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB) override
//     {
// // 	  GetTimeTable()->Start(Info());

//         BaseType::GetPreconditioner()->Initialize(rA,rX,rB);

//         bool is_solved = true;
//         VectorType x(TDenseSpaceType::Size1(rX));
//         VectorType b(TDenseSpaceType::Size1(rB));
//         for(unsigned int i = 0 ; i < TDenseSpaceType::Size2(rX) ; i++)
//         {
//             TDenseSpaceType::GetColumn(i,rX, x);
//             TDenseSpaceType::GetColumn(i,rB, b);

//             BaseType::GetPreconditioner()->ApplyInverseRight(x);
//             BaseType::GetPreconditioner()->ApplyLeft(b);

//             is_solved &= IterativeSolve(rA,x,b);

//             BaseType::GetPreconditioner()->Finalize(x);
//         }

// // 	  GetTimeTable()->Stop(Info());

//         return is_solved;
//     }

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
        const int size = rX.size();

        BaseType::mIterationsNumber = 0;

        VectorType r(size);
        this->PreconditionedMult(rpLinearOperator, rX, r);
        r *= -1.0;
        r += rB;

        BaseType::mBNorm = rB.Norm();

        VectorType p(r);
        VectorType q(size);
        q.SetValue(0.0);

        double roh0 = r.Dot(r);
        double roh1 = roh0;
        double beta = 0.0;

        if(std::abs(roh0) < std::numeric_limits<double>::epsilon()) {
            return false;
        }

        do
        {
            this->PreconditionedMult(rpLinearOperator, p, q);

            double pq = p.Dot(q);

            if(std::abs(pq) < std::numeric_limits<double>::epsilon()) {
                break;
            }

            double alpha = roh0 / pq;

            rX.Add(alpha, p);
            r.Add(-alpha, q);

            roh1 = r.Dot(r);

            beta = (roh1 / roh0);
            p *= beta;
            p.Add(1.0, r);

            roh0 = roh1;

            BaseType::mResidualNorm = std::sqrt(roh1);
            BaseType::mIterationsNumber++;
        }
        while(BaseType::IterationNeeded() && (std::abs(roh0) > std::numeric_limits<double>::epsilon()));

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



