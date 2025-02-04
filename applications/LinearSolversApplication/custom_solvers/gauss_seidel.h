//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

#pragma once

// Project includes
#include "includes/smart_pointers.h"
#include "linear_solvers/linear_solver.h" // LinearSolver
#include "factories/linear_solver_factory.h" // LinearSolverFactory
#include "includes/code_location.h" // KRATOS_CODE_LOCATION
#include "includes/define.h" // KRATOS_CLASS_POINTER_DEFINITION

// System includes
#include <memory> // std::unique_ptr


namespace Kratos {


template <class TSparseSpace, class TDenseSpace>
class KRATOS_API(LINEARSOLVERS_APPLICATION) GaussSeidelSmoother
    : public LinearSolver<TSparseSpace,TDenseSpace,Reorderer<TSparseSpace,TDenseSpace>>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(GaussSeidelSmoother);

    using Base = LinearSolver<TSparseSpace,TDenseSpace,Reorderer<TSparseSpace,TDenseSpace>>;

    using typename Base::SparseMatrixType;

    using typename Base::VectorType;

    /// @copydoc Base::Base
    GaussSeidelSmoother();

    GaussSeidelSmoother(Parameters Settings);

    /// @details Required by PIMPL.
    ~GaussSeidelSmoother() override;

    /// @copydoc Base::Initialize
    void Initialize(SparseMatrixType& rLhs,
                    VectorType& rSolution,
                    VectorType& rRhs) override;

    /// @copydoc Base::InitializeSolutionStep
    void InitializeSolutionStep(SparseMatrixType& rLhs,
                                VectorType& rSolution,
                                VectorType& rRhs) override;

    bool Solve(SparseMatrixType& rLhs,
               VectorType& rSolution,
               VectorType& rRhs) override;

    /// @copydoc Base::FinalizeSolutionStep
    void FinalizeSolutionStep(SparseMatrixType& rLhs,
                              VectorType& rSolution,
                              VectorType& rRhs) override;

    Parameters GetDefaultParameters() const;

    /// @name Unsupported Interface
    /// @{

    /// @copydoc Base::SetTolerance
    void SetTolerance(double NewTolerance) override
    {KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented\n";}

    /// @}

private:
    struct Impl;
    std::unique_ptr<Impl> mpImpl;
}; // class GaussSeidelSmoother


template <class TSparseSpace, class TDenseSpace>
class KRATOS_API(LINEARSOLVERS_APPLICATION) GaussSeidelSmootherFactory
    : public LinearSolverFactory<TSparseSpace,TDenseSpace>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(GaussSeidelSmootherFactory);

protected:
    typename LinearSolver<TSparseSpace,TDenseSpace>::Pointer CreateSolver(Parameters Settings) const override;
}; // class GaussSeidelSmootherFactory


} // namespace Kratos
