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

    using typename Base::Status;

    NoOpConstraintAssembler()
        : NoOpConstraintAssembler(Parameters())
    {}

    NoOpConstraintAssembler(Parameters Settings)
        : NoOpConstraintAssembler(Settings, "unnamed")
    {}

    NoOpConstraintAssembler(Parameters, std::string&& rInstanceName) noexcept
        : Base(ConstraintImposition::None, std::move(rInstanceName))
    {}

    Status FinalizeSolutionStep(typename TSparse::MatrixType& rLhs,
                                typename TSparse::VectorType& rSolution,
                                typename TSparse::VectorType& rRhs,
                                const std::size_t iIteration) override
    {
        return Status {/*finished=*/true, /*converged=*/true};
    }
}; // class NoOpConstraintAssembler


} // namespace Kratos

