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
class MasterSlaveConstraintAssembler : public ConstraintAssembler<TSparse,TDense>
{
private:
    using Base = ConstraintAssembler<TSparse,TDense>;

public:
    MasterSlaveConstraintAssembler() noexcept;

    MasterSlaveConstraintAssembler(Parameters Settings);

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
                  const typename Base::DofSet& rDofSet,
                  const bool AssembleLhs,
                  const bool AssembleRhs) override;

    /// @copydoc Base::Initialize
    void Initialize(typename TSparse::MatrixType& rLhs,
                    typename TSparse::VectorType& rRhs) override;

    /// @copydoc Base::FinalizeSolutionStep
    typename Base::Status FinalizeSolutionStep(typename TSparse::MatrixType& rLhs,
                                               typename TSparse::VectorType& rSolution,
                                               typename TSparse::VectorType& rRhs,
                                               const std::size_t iIteration) override;

    /// @copydoc Base::Finalize
    void Finalize(typename TSparse::MatrixType& rLhs,
                  typename TSparse::VectorType& rSolution,
                  typename TSparse::VectorType& rRhs,
                  typename Base::DofSet& rDofSet) override;

    /// @copydoc Base::Clear
    void Clear() override;

    static Parameters GetDefaultParameters();

private:
    std::vector<std::size_t> mSlaveIds;

    std::vector<std::size_t> mMasterIds;

    std::unordered_set<std::size_t> mInactiveSlaveIds;

    DiagonalScaling mDiagonalScaling;
}; // class MasterSlaveConstraintAssembler


} // namespace Kratos
