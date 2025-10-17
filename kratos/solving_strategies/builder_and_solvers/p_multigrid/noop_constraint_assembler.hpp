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


/// @brief Dummy constraint assembler class that does not impose constraints.
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
    {
    }

    NoOpConstraintAssembler(Parameters Settings, std::string&& rInstanceName) noexcept
        : Base(ConstraintImposition::None, std::move(rInstanceName)),
          mVerbosity(Settings.Has("verbosity") ? Settings["verbosity"].Get<int>() : 1)
    {}

    bool FinalizeSolutionStep(typename TSparse::MatrixType& rLhs,
                              typename TSparse::VectorType& rSolution,
                              typename TSparse::VectorType& rRhs,
                              PMGStatusStream::Report& rReport,
                              PMGStatusStream& rStream) override
    {
        rReport.maybe_constraint_residual = 0;
        rReport.constraints_converged = true;
        rStream.Submit(rReport.Tag(2), mVerbosity);
        return true;
    }

private:
    int mVerbosity;
}; // class NoOpConstraintAssembler


} // namespace Kratos

