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
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "future/containers/linear_system_tags.h"
#include "future/linear_operators/linear_operator.h"
#include "future/linear_solvers/linear_solver.h"
#include "future/preconditioners/preconditioner.h"

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
 * @class IterativeSolver
 * @ingroup KratosCore
 * @brief Base class for iterative linear solvers.
 * @details This class serves as a base for all iterative linear solvers.
 * It defines the common interface and data members for iterative solvers,
 * such as the tolerance, maximum number of iterations, and preconditioner.
 * @tparam TLinearAlgebra Type of the linear algebra used.
 */
template<class TLinearAlgebra>
class IterativeSolver : public Future::LinearSolver<TLinearAlgebra>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of IterativeSolver
    KRATOS_CLASS_POINTER_DEFINITION(IterativeSolver);

    /// The base class definition
    using BaseType = Future::LinearSolver<TLinearAlgebra>;

    /// Type definition for index
    using IndexType = typename TLinearAlgebra::IndexType;

    /// Vector type definition from linear algebra template parameter
    using VectorType = typename TLinearAlgebra::VectorType;

    /// Linear operator type definition
    using LinearOperatorType = LinearOperator<TLinearAlgebra>;

    /// Preconditioner pointer type definition
    using PreconditionerPointerType = typename Preconditioner<TLinearAlgebra>::Pointer; //TODO: maybe it is a good idea to make this unique_ptr

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    IterativeSolver()
        : BaseType()
    {}

    /// Copy constructor.
    IterativeSolver(const IterativeSolver& Other) = delete;

    /// Destructor.
    ~IterativeSolver() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    IterativeSolver& operator=(const IterativeSolver& Other) = delete;

    ///@}
    ///@name Operations
    ///@{

    void Initialize(LinearSystem<TLinearAlgebra>& rLinearSystem) override
    {
        // Call the preconditioner initialize
        const auto lhs_tag = LinearSystemTags::SparseMatrixTagFromString(BaseType::mLhsTagString);
        const auto& rp_lhs_lin_op = rLinearSystem.pGetLinearOperator(lhs_tag);
        this->GetPreconditioner()->Initialize(rp_lhs_lin_op);
    }

    void InitializeSolutionStep(LinearSystem<TLinearAlgebra>& rLinearSystem) override
    {
        // Call the preconditioner initialize solution step
        const auto lhs_tag = LinearSystemTags::SparseMatrixTagFromString(BaseType::mLhsTagString);
        const auto& rp_lhs_lin_op = rLinearSystem.pGetLinearOperator(lhs_tag);
        this->GetPreconditioner()->InitializeSolutionStep(rp_lhs_lin_op);
    }

    bool PerformSolutionStep(LinearSystem<TLinearAlgebra>& rLinearSystem) override
    {
        bool is_solved;
        if (!this->mMultipleSolve) {
            // Get sparse matrix and dense vector tags from strings
            const auto dx_tag = LinearSystemTags::DenseVectorTagFromString(BaseType::mDxTagString);
            const auto rhs_tag = LinearSystemTags::DenseVectorTagFromString(BaseType::mRhsTagString);
            const auto lhs_tag = LinearSystemTags::SparseMatrixTagFromString(BaseType::mLhsTagString);

            // Check if the linear system is consistent
            if(!rLinearSystem.IsConsistent(lhs_tag, rhs_tag, dx_tag)) {
                KRATOS_WARNING("IterativeSolver") << "Linear system is not consistent. PerformSolutionStep cannot be performed." << std::endl;
                return false;
            }

            // Get arrays from linear system container
            auto& r_dx = *(rLinearSystem.pGetVector(dx_tag));
            auto& r_rhs = *(rLinearSystem.pGetVector(rhs_tag));
            const auto& rp_lhs_lin_op = rLinearSystem.pGetLinearOperator(lhs_tag);

            // Call CG implementation
            is_solved = this->IterativeSolve(rp_lhs_lin_op, r_rhs, r_dx);
            KRATOS_WARNING_IF("IterativeSolver", !is_solved) << "Non converged linear solution. ["<< GetResidualNorm() / mBNorm << " > "<<  GetTolerance() << "]" << std::endl;

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
                KRATOS_WARNING("IterativeSolver") << "Linear system is not consistent. PerformSolutionStep cannot be performed." << std::endl;
                return false;
            }

            // Perform the columnwise solve
            is_solved = true;
            VectorType aux_dx(r_rhs.size1()); // Auxiliary update vector
            VectorType aux_rhs(r_rhs.size1()); // Auxiliary residual vector
            const std::size_t n_problems = r_rhs.size2();

            KRATOS_WATCH(r_rhs.size1())
            KRATOS_WATCH(r_rhs.size2())
            for (std::size_t i = 0; i < n_problems; ++i) {
                // Get current residual column
                IndexPartition<IndexType>(aux_rhs.size()).for_each([i, &r_rhs, &aux_rhs](const IndexType Index) {
                    aux_rhs[Index] = r_rhs(Index, i);
                });

                // Call CG implementation
                const bool aux_is_solved = this->IterativeSolve(rp_lhs_lin_op, aux_rhs, aux_dx);
                KRATOS_WARNING_IF("IterativeSolver", !aux_is_solved) << "Non converged linear solution for column " << i << ". ["<< GetResidualNorm() / mBNorm << " > "<<  GetTolerance() << "]" << std::endl;
                is_solved &= aux_is_solved;

                // Set current solution update
                IndexPartition<IndexType>(aux_dx.size()).for_each([i, &r_dx, &aux_dx](const IndexType Index) {
                    r_dx(Index, i) = aux_dx[Index];
                });
            }
        }

        return is_solved;
    }

    void FinalizeSolutionStep(LinearSystem<TLinearAlgebra>& rLinearSystem) override
    {
        // Call the preconditioner finalize solution step
        const auto lhs_tag = LinearSystemTags::SparseMatrixTagFromString(BaseType::mLhsTagString);
        const auto& rp_lhs_lin_op = rLinearSystem.pGetLinearOperator(lhs_tag);
        this->GetPreconditioner()->FinalizeSolutionStep(rp_lhs_lin_op);
    }

    void Clear() override
    {
        this->GetPreconditioner()->Clear();
    }

    ///@}
    ///@name Access
    ///@{

    virtual PreconditionerPointerType GetPreconditioner()
    {
        KRATOS_ERROR_IF_NOT(mpPreconditioner) << "Preconditioner not set" << std::endl;
        return mpPreconditioner;
    }

    virtual const PreconditionerPointerType GetPreconditioner() const
    {
        KRATOS_ERROR_IF_NOT(mpPreconditioner) << "Preconditioner not set" << std::endl;
        return mpPreconditioner;
    }

    virtual void SetPreconditioner(PreconditionerPointerType pNewPreconditioner)
    {
        mpPreconditioner = pNewPreconditioner;
    }

    virtual void SetMaxIterationsNumber(unsigned int NewMaxIterationsNumber)
    {
        mMaxIterationsNumber = NewMaxIterationsNumber;
    }

    virtual IndexType GetMaxIterationsNumber()
    {
        return mMaxIterationsNumber;
    }

    virtual void SetIterationsNumber(unsigned int NewIterationNumber)
    {
        mIterationsNumber = NewIterationNumber;
    }

    IndexType GetIterationsNumber() override
    {
        return mIterationsNumber;
    }

    void SetTolerance(double NewTolerance) override
    {
        mTolerance = NewTolerance;
    }

    double GetTolerance() override
    {
        return mTolerance;
    }

    virtual void SetResidualNorm(double NewResidualNorm)
    {
        if (mIterationsNumber == 1)
            mFirstResidualNorm = NewResidualNorm;
        mResidualNorm = NewResidualNorm;
    }

    virtual double GetResidualNorm()
    {
        return mResidualNorm;
    }

    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters( R"({
            "solver_type" : "iterative_solver",
            "tolerance" : 1e-6,
            "max_iteration" : 100,
            "preconditioner_type" : "none"
        })");
        default_parameters.AddMissingParameters(BaseType::GetDefaultParameters());

        return default_parameters;
    }

    ///@}
    ///@name Inquiry
    ///@{

    virtual bool IterationNeeded()
    {
        return (mIterationsNumber < mMaxIterationsNumber) && (mResidualNorm > mTolerance * mBNorm);
    }

    virtual bool IsConverged()
    {
        return (mResidualNorm <= mTolerance * mBNorm);
    }

    bool RequiresAdditionalData() const override
    {
        if (GetPreconditioner()->RequiresAdditionalData())
            return true;
        else
            return false;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Iterative solver with " << GetPreconditioner()->Info();
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
        if (mBNorm < std::numeric_limits<double>::epsilon()) {
            if (mResidualNorm > std::numeric_limits<double>::epsilon()) {
                rOStream << "    Residual ratio : infinite" << std::endl;
            } else {
                rOStream << "    Residual ratio : 0.0" << std::endl;
            }
        } else {
            rOStream << "    Initial Residual ratio : " << mBNorm << std::endl;
            rOStream << "    Final Residual ratio : " << mResidualNorm << std::endl;
            rOStream << "    Residual ratio : " << mResidualNorm / mBNorm << std::endl;
            rOStream << "    Slope : " << (mResidualNorm - mBNorm) / mIterationsNumber << std::endl;
        }

        rOStream << "    Tolerance : " << mTolerance << std::endl;
        rOStream << "    Number of iterations : " << mIterationsNumber << std::endl;
        rOStream << "    Maximum number of iterations : " << mMaxIterationsNumber << std::endl;
        KRATOS_WARNING_IF("IterativeSolver", mMaxIterationsNumber == mIterationsNumber) << "Iterative solver non converged! " << std::endl;
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

    double mBNorm = 0.0;

    double mResidualNorm = 0.0;

    double mFirstResidualNorm = 0.0;

    IndexType mIterationsNumber = 0;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    virtual bool IterativeSolve(
        const typename LinearOperatorType::UniquePointer& rpLinearOperator,
        const VectorType& rB,
        VectorType& rX)
    {
        KRATOS_ERROR << "Calling IterativeSolve in IterativeSolver base class. Implement it in derived class." << std::endl;
        return false;
    }

    void AssignSettings(const Parameters& Settings) override
    {
        // Assign base class settings
        BaseType::AssignSettings(Settings);

        // Assign input settings to member variables
        mTolerance = Settings["tolerance"].GetDouble();
        mMaxIterationsNumber = Settings["max_iteration"].GetInt();
    }

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

    double mTolerance = 0.0;

    IndexType mMaxIterationsNumber = 0;

    PreconditionerPointerType mpPreconditioner = nullptr;

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

}; // Class IterativeSolver

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
    IterativeSolver<TLinearAlgebra>& rThis)
{
    return IStream;
}

/// output stream function
template<class TLinearAlgebra>
inline std::ostream& operator << (
    std::ostream& OStream,
    const IterativeSolver<TLinearAlgebra>& rThis)
{
    rThis.PrintInfo(OStream);
    OStream << std::endl;
    rThis.PrintData(OStream);

    return OStream;
}

///@}

}  // namespace Kratos::Future.
