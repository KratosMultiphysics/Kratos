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
private:
    using Base = ConstraintAssembler<TSparse,TDense>;

public:
    AugmentedLagrangeConstraintAssembler() noexcept;

    AugmentedLagrangeConstraintAssembler(Parameters Settings);

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
                  const typename Base::DofSet& rDofSet) override;

    /// @copydoc Base::Initialize
    void Initialize(const typename Base::ConstraintArray& rConstraints,
                    const ProcessInfo& rProcessInfo,
                    typename TSparse::MatrixType& rLhs,
                    typename TSparse::VectorType& rRhs,
                    typename Base::DofSet& rDofSet) override;

    /// @copydoc Base::FinalizeSolutionStep
    typename Base::Status FinalizeSolutionStep(const typename Base::ConstraintArray& rConstraints,
                                               const ProcessInfo& rProcessInfo,
                                               typename TSparse::MatrixType& rLhs,
                                               typename TSparse::VectorType& rSolution,
                                               typename TSparse::VectorType& rRhs,
                                               const typename Base::DofSet& rDofSet,
                                               const std::size_t iIteration) override;

    /// @copydoc Base::Clear
    void Clear() override;

    static Parameters GetDefaultParameters();

    typename TSparse::DataType GetPenaltyFactor() const
    {
        return this->GetValue(AugmentedLagrangeConstraintAssembler::GetAlgorithmicParametersVariable())[0];
    }

    typename TSparse::DataType GetInitialLagrangeMultiplier() const
    {
        return this->GetValue(AugmentedLagrangeConstraintAssembler::GetAlgorithmicParametersVariable())[1];
    }

    typename TSparse::DataType GetTolerance() const
    {
        return this->GetValue(AugmentedLagrangeConstraintAssembler::GetAlgorithmicParametersVariable())[2];
    }

    static const Variable<Vector>& GetAlgorithmicParametersVariable() noexcept
    {
        return SHAPE_FUNCTIONS_VECTOR;
    }

private:
    /// @brief A map associating slave IDs with constraint indices.
    std::unordered_map<std::size_t,std::size_t> mSlaveToConstraintMap;

    typename TSparse::MatrixType mTransposeRelationMatrix;

    int mVerbosity;
}; // class AugmentedLagrangeConstraintAssembler


} // namespace Kratos
