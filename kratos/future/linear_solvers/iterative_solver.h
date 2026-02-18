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
#include "future/linear_solvers/preconditioner.h"

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

/// Base class for all the iterative solvers in Kratos.
/** This class define the general interface for the iterative solvers in Kratos.
    iterative solver is a template class with this parameter:
    - TSparseSpaceType which specify type
      of the unknowns, coefficients, sparse matrix, vector of
  unknowns, right hand side vector and their respective operators.
    - TDenseMatrixType which specify type of the
      matrices used as temporary matrices or multi solve unknowns and
  right hand sides and their operators.
    - TPreconditionerType  which specify type of the preconditioner to be used.
    - TStopCriteriaType for specifying type of the object which control the stop criteria for iteration loop.
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

    /// Type definition for data
    using DataType = typename TLinearAlgebra::DataType;

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
    IterativeSolver(Parameters Settings = Parameters(R"({})"))
        : BaseType(Settings)
    {
        // Validate and assign default parameters
        Settings.ValidateAndAssignDefaults(this->GetDefaultParameters());

        // Assign input settings to member variables
        mTolerance = Settings["tolerance"].GetDouble();
        mMaxIterationsNumber = Settings["max_iteration"].GetInt();
        mpPreconditioner = Kratos::make_shared<Preconditioner<TLinearAlgebra>>(); //TODO: implement preconditioner by leveraging the registry

        // Assign the linear system tags to be used
        this->mDxTagString = Settings["dx_tag"].GetString();
        this->mRhsTagString = Settings["rhs_tag"].GetString();
        this->mLhsTagString = Settings["lhs_tag"].GetString();
    }

    /// Constructor with preconditioner.
    IterativeSolver(
        Parameters Settings,
        PreconditionerPointerType pPreconditioner)
        : BaseType(Settings)
    {
        Settings.ValidateAndAssignDefaults(this->GetDefaultParameters());
        SetPreconditioner(pPreconditioner);
    }

    // IterativeSolver(Parameters settings,
    //                 typename TPreconditionerType::Pointer pNewPreconditioner = Kratos::make_shared<TPreconditionerType>()
    //                ):
    //     mResidualNorm(0),
    //     mIterationsNumber(0),
    //     mBNorm(0),
    //     mpPreconditioner(pNewPreconditioner)
    // {
    //     KRATOS_TRY

    //     Parameters default_parameters( R"(
    //     {
    //     "solver_type": "IterativeSolver",
    //     "tolerance" : 1.0e-6,
    //     "max_iteration" : 200,
    //     "preconditioner_type": "none",
    //     "scaling":false
    //     }  )" );

    //     //now validate agains defaults -- this also ensures no type mismatch
    //     settings.ValidateAndAssignDefaults(default_parameters);

    //     this->SetTolerance( settings["tolerance"].GetDouble() );
    //     this->SetMaxIterationsNumber( settings["max_iteration"].GetInt() );


    //     KRATOS_CATCH("")
    // }


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
        // Get sparse matrix and dense vector tags from strings
        const auto dx_tag = LinearSystemTags::DenseVectorTagFromString(BaseType::mDxTagString);
        const auto rhs_tag = LinearSystemTags::DenseVectorTagFromString(BaseType::mRhsTagString);
        const auto lhs_tag = LinearSystemTags::SparseMatrixTagFromString(BaseType::mLhsTagString);

        // Call the preconditioner initialize
        const auto& rp_lhs_lin_op = rLinearSystem.pGetLinearOperator(lhs_tag);
        auto p_rhs = rLinearSystem.pGetVector(rhs_tag);
        auto p_dx = rLinearSystem.pGetVector(dx_tag);
        this->GetPreconditioner()->Initialize(rp_lhs_lin_op, p_rhs, p_dx);
    }

    void InitializeSolutionStep(LinearSystem<TLinearAlgebra>& rLinearSystem) override
    {
        // Get sparse matrix and dense vector tags from strings
        const auto dx_tag = LinearSystemTags::DenseVectorTagFromString(BaseType::mDxTagString);
        const auto rhs_tag = LinearSystemTags::DenseVectorTagFromString(BaseType::mRhsTagString);
        const auto lhs_tag = LinearSystemTags::SparseMatrixTagFromString(BaseType::mLhsTagString);

        // Call the preconditioner initialize
        const auto& rp_lhs_lin_op = rLinearSystem.pGetLinearOperator(lhs_tag);
        auto p_rhs = rLinearSystem.pGetVector(rhs_tag);
        auto p_dx = rLinearSystem.pGetVector(dx_tag);
        this->GetPreconditioner()->InitializeSolutionStep(rp_lhs_lin_op, p_rhs, p_dx);
    }

    void FinalizeSolutionStep(LinearSystem<TLinearAlgebra>& rLinearSystem) override
    {
        // Get sparse matrix and dense vector tags from strings
        const auto dx_tag = LinearSystemTags::DenseVectorTagFromString(BaseType::mDxTagString);
        const auto rhs_tag = LinearSystemTags::DenseVectorTagFromString(BaseType::mRhsTagString);
        const auto lhs_tag = LinearSystemTags::SparseMatrixTagFromString(BaseType::mLhsTagString);

        // Call the preconditioner finalize
        const auto& rp_lhs_lin_op = rLinearSystem.pGetLinearOperator(lhs_tag);
        auto p_rhs = rLinearSystem.pGetVector(rhs_tag);
        auto p_dx = rLinearSystem.pGetVector(dx_tag);
        this->GetPreconditioner()->FinalizeSolutionStep(rp_lhs_lin_op, p_rhs, p_dx);
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
            "dx_tag" : "Dx",
            "rhs_tag" : "RHS",
            "lhs_tag" : "LHS",
            "tolerance" : 1e-6,
            "max_iteration" : 100,
            "multiple_solve" : false
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
        if (mBNorm == 0.00)
            if (mResidualNorm != 0.00)
                rOStream << "    Residual ratio : infinite" << std::endl;
            else
                rOStream << "    Residual ratio : 0" << std::endl;
        else
        {
            rOStream << "    Initial Residual ratio : " << mBNorm << std::endl;
            rOStream << "    Final Residual ratio : " << mResidualNorm << std::endl;
            rOStream << "    Residual ratio : " << mResidualNorm / mBNorm << std::endl;
            rOStream << "    Slope : " << (mResidualNorm - mBNorm) / mIterationsNumber << std::endl;
        }

        rOStream << "    Tolerance : " << mTolerance << std::endl;
        rOStream << "    Number of iterations : " << mIterationsNumber << std::endl;
        rOStream << "    Maximum number of iterations : " << mMaxIterationsNumber;
        if (mMaxIterationsNumber == mIterationsNumber)
            rOStream << std::endl << "!!!!!!!!!!!! ITERATIVE SOLVER NON CONVERGED !!!!!!!!!!!!" << mMaxIterationsNumber;
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

    void PreconditionedMult(
        const typename LinearOperatorType::UniquePointer& rpLinearOperator,
        const VectorType& rX,
        VectorType& rY)
    {
        GetPreconditioner()->Mult(rpLinearOperator, rX, rY);
    }

    void PreconditionedTransposeMult(
        const typename LinearOperatorType::UniquePointer& rpLinearOperator,
        const VectorType& rX,
        VectorType& rY)
    {
        GetPreconditioner()->TransposeMult(rpLinearOperator, rX, rY);
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
