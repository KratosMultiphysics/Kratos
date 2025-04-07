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
class AugmentedLagrangeConstraintAssembler : public ConstraintAssembler<TSparse,TDense>
{
public:
    using Base = ConstraintAssembler<TSparse,TDense>;

    AugmentedLagrangeConstraintAssembler() noexcept;

    AugmentedLagrangeConstraintAssembler(Parameters Settings);

    AugmentedLagrangeConstraintAssembler(Parameters Settings,
                                         std::string&& rInstanceName);

    ~AugmentedLagrangeConstraintAssembler();

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

    /// @copydoc Base::InitializeSolutionStep
    void InitializeSolutionStep(typename TSparse::MatrixType& rLhs,
                                typename TSparse::VectorType& rSolution,
                                typename TSparse::VectorType& rRhs) override;

    /// @copydoc Base::FinalizeSolutionStep
    bool FinalizeSolutionStep(typename TSparse::MatrixType& rLhs,
                              typename TSparse::VectorType& rSolution,
                              typename TSparse::VectorType& rRhs,
                              PMGStatusStream::Report& rReport) override;

    /// @copydoc Base::Finalize
    void Finalize(typename TSparse::MatrixType& rLhs,
                  typename TSparse::VectorType& rSolution,
                  typename TSparse::VectorType& rRhs,
                  typename Base::DofSet& rDofSet) override;

    /// @copydoc Base::Clear
    void Clear() override;

    static Parameters GetDefaultParameters();

    typename TSparse::DataType GetPenaltyFactor() const;

    typename TSparse::DataType GetInitialLagrangeMultiplier() const
    {
        return this->GetValue(AugmentedLagrangeConstraintAssembler::GetAlgorithmicParametersVariable())[0];
    }

    typename TSparse::DataType GetTolerance() const
    {
        return this->GetValue(AugmentedLagrangeConstraintAssembler::GetAlgorithmicParametersVariable())[1];
    }

    static const Variable<Vector>& GetAlgorithmicParametersVariable() noexcept
    {
        return SHAPE_FUNCTIONS_VECTOR;
    }

private:
    typename TSparse::MatrixType& GetTransposeRelationMatrix();

    struct Impl;
    std::unique_ptr<Impl> mpImpl;
}; // class AugmentedLagrangeConstraintAssembler


} // namespace Kratos
