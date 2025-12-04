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
class MasterSlaveConstraintAssembler final : public ConstraintAssembler<TSparse,TDense>
{
public:
    using Base = ConstraintAssembler<TSparse,TDense>;

    MasterSlaveConstraintAssembler() noexcept;

    MasterSlaveConstraintAssembler(Parameters Settings);

    MasterSlaveConstraintAssembler(Parameters Settings,
                                   std::string&& rInstanceName);

    ~MasterSlaveConstraintAssembler();

    /// @copydoc Base::AllocateConstraints
    void AllocateConstraints(PointerVectorSet<MasterSlaveConstraint,IndexedObject>::const_iterator itConstraintBegin,
                             PointerVectorSet<MasterSlaveConstraint,IndexedObject>::const_iterator itConstraintEnd,
                             const ProcessInfo& rProcessInfo,
                             typename Base::DofSet::const_iterator itDofBegin,
                             typename Base::DofSet::const_iterator itDofEnd) override;

    /// @copydoc Base::Assemble
    void Assemble(const typename Base::ConstraintArray& rConstraints,
                  const ProcessInfo& rProcessInfo,
                  typename Base::DofSet::const_iterator itDofBegin,
                  typename Base::DofSet::const_iterator itDofEnd,
                  const bool AssembleLhs,
                  const bool AssembleRhs) override;

    /// @copydoc Base::Initialize
    void Initialize(typename TSparse::MatrixType& rLhs,
                    typename TSparse::VectorType& rSolution,
                    typename TSparse::VectorType& rRhs,
                    typename Base::DofSet& rDofs) override;

    /// @copydoc Base::FinalizeSolutionStep
    bool FinalizeConstraintIteration(typename TSparse::MatrixType& rLhs,
                                     typename TSparse::VectorType& rSolution,
                                     typename TSparse::VectorType& rRhs,
                                     typename Base::DofSet::iterator itDofBegin,
                                     typename Base::DofSet::iterator itDofEnd,
                                     PMGStatusStream::Report& rReport,
                                     PMGStatusStream& rStream) override;

    /// @copydoc Base::Finalize
    void Finalize(typename TSparse::MatrixType& rLhs,
                  typename TSparse::VectorType& rSolution,
                  typename TSparse::VectorType& rRhs,
                  typename Base::DofSet& rDofSet) override;

    /// @copydoc Base::ComputeDependentResidual
    void ComputeIndependentResidual(typename TSparse::VectorType& rResidual) const override;

    /// @copydoc Base::ComputeDependentResidual
    void ComputeDependentResidual(typename TSparse::VectorType& rResidual) const override;

    /// @copydoc Base::ComputeIndependentSolution
    void ComputeIndependentSolution(typename TSparse::VectorType& rSolution) const override;

    /// @copydoc Base::ComputeDependentSolution
    void ComputeDependentSolution(typename TSparse::VectorType& rSolution) const override;

    /// @copydoc Base::GetDependentDofs
    const typename Base::DofSet& GetDependentDofs(const typename Base::DofSet& rIndependentDofSet) const noexcept override;

    /// @copydoc Base::Clear
    void Clear() override;

    static Parameters GetDefaultParameters();

private:
    struct Impl;
    std::unique_ptr<Impl> mpImpl;
}; // class MasterSlaveConstraintAssembler


} // namespace Kratos
