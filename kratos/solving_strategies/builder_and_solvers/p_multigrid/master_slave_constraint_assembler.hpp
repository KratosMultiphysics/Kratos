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
#include "solving_strategies/builder_and_solvers/p_multigrid/diagonal_scaling.hpp" // ParseDiagonalScaling, GetDiagonalScaleFactor
#include "includes/kratos_parameters.h" // Parameters


namespace Kratos {


template <class TSparse, class TDense>
class MasterSlaveConstraintAssembler final : public ConstraintAssembler<TSparse,TDense>
{
public:
    using Base = ConstraintAssembler<TSparse,TDense>;

    MasterSlaveConstraintAssembler() noexcept;

    MasterSlaveConstraintAssembler(Parameters Settings);

    MasterSlaveConstraintAssembler(Parameters Settings,
                                   std::string&& rInstanceName);

    ~MasterSlaveConstraintAssembler();

    /// @copydoc Base::Allocate
    void Allocate(const typename Base::ConstraintArray& rConstraints,
                  const ProcessInfo& rProcessInfo,
                  typename TSparse::MatrixType& rLhs,
                  typename TSparse::VectorType& rSolution,
                  typename TSparse::VectorType& rRhs,
                  typename Base::DofSet& rDofSet) override;

    /// @copydoc Base::Assemble
    void Assemble(const typename Base::ConstraintArray& rConstraints,
                  const ProcessInfo& rProcessInfo,
                  typename Base::DofSet& rDofSet,
                  const bool AssembleLhs,
                  const bool AssembleRhs) override;

    /// @copydoc Base::Initialize
    void Initialize(typename TSparse::MatrixType& rLhs,
                    typename TSparse::VectorType& rRhs,
                    typename Base::DofSet::iterator itDofBegin,
                    typename Base::DofSet::iterator itDofEnd) override;

    /// @copydoc Base::FinalizeSolutionStep
    bool FinalizeSolutionStep(typename TSparse::MatrixType& rLhs,
                              typename TSparse::VectorType& rSolution,
                              typename TSparse::VectorType& rRhs,
                              PMGStatusStream::Report& rReport) override;

    /// @copydoc Base::Apply
    void Apply(typename TSparse::VectorType& rSolution) const override;

    /// @copydoc Base::Clear
    void Clear() override;

    static Parameters GetDefaultParameters();

private:
    struct Impl;
    std::unique_ptr<Impl> mpImpl;
}; // class MasterSlaveConstraintAssembler


} // namespace Kratos
