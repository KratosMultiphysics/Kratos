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

// System includes
#include <optional> // std::optional


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
    typename TSparse::MatrixType& GetTransposeRelationMatrix();

    /// @brief A map associating slave IDs with constraint indices.
    std::unordered_map<std::size_t,std::size_t> mSlaveToConstraintMap;

    std::optional<typename TSparse::MatrixType> mMaybeTransposeRelationMatrix;

    int mVerbosity;
}; // class AugmentedLagrangeConstraintAssembler


} // namespace Kratos
