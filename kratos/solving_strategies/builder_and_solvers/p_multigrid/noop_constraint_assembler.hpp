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

// Project Includes
#include "solving_strategies/builder_and_solvers/p_multigrid/constraint_assembler.hpp" // ConstraintAssembler
#include "includes/kratos_parameters.h" // Parameters


namespace Kratos {


template <class TSparse, class TDense>
class NoOpConstraintAssembler final : public ConstraintAssembler<TSparse,TDense>
{
public:
    using Base = ConstraintAssembler<TSparse,TDense>;

    NoOpConstraintAssembler()
        : NoOpConstraintAssembler(Parameters())
    {}

    NoOpConstraintAssembler(Parameters Settings)
        : NoOpConstraintAssembler(Settings, "unnamed")
    {}

    NoOpConstraintAssembler(Parameters, std::string&& rInstanceName) noexcept
        : Base(ConstraintImposition::None, std::move(rInstanceName))
    {}

    bool FinalizeSolutionStep(typename TSparse::MatrixType& rLhs,
                                typename TSparse::VectorType& rSolution,
                                typename TSparse::VectorType& rRhs,
                                PMGStatusStream::Report& rReport) override
    {
        rReport.maybe_constraint_residual = 0;
        rReport.constraints_converged = true;
        return true;
    }
}; // class NoOpConstraintAssembler


} // namespace Kratos

