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


/// @brief Class implementing constraint imposition with the augmented Lagrange multiplier method.
/// @details In a nutshell, the augmented Lagrange method is a combination of the penalty method
///          and the method of Lagrange multipliers. The solution potentially involves multiple
///          iterations, depending on the choice of penalty factor as well as the properties of
///          the system to be solved.
///
///          Default settings:
///          @code
///          {
///             "method" : "augmented_lagrange",
///             "penalty_factor" : "norm",
///             "initial_lagrange_multiplier" : 0.0,
///             "tolerance" : 1e-6,
///             "max_iterations" : 1e1,
///             "verbosity" : 1
///          }
///          @endcode
///          - @p "method" name of the constrained assembler to refer to in @ref ConstraintAssemblerFactory.
///          - @p "penalty_factor" Coefficient to scale the constraint equations by. It can either be a
///                                numeric literal or a string representing a function of the left hand side
///                                matrix' diagonal. See @ref Scaling for more information.
///             - @p "max" Infinity norm (abs max) of the matrix' diagonal.
///             - @p "norm" 2-norm of the matrix' diagonal.
///         - @p "initial_lagrange_multiplier" Value to initialize lagrange multipliers with.
///         - @p "tolerance" Target absolute norm of the constraint equations (without penalty scaling).
///         - @p "max_iterations" Stop iterating if the desired tolerance was not reached after this many iterations.
///         - @p "verbosity" Level of verbosity. Every level includes lower verbosity levels and
///                          adds new events to report:
///             - @p 0 No messages, not even in case of failure.
///             - @p 1 Warnings and failure messages.
///             - @p 2 Aggregated status reports.
///             - @p 3 Per-iteration status reports.
///             - @p 4 Output system matrices and vectors.
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

    /// @internal
    typename TSparse::DataType GetPenaltyFactor() const;

    /// @internal
    typename TSparse::DataType GetInitialLagrangeMultiplier() const
    {
        return this->GetValue(AugmentedLagrangeConstraintAssembler::GetAlgorithmicParametersVariable())[0];
    }

    /// @internal
    typename TSparse::DataType GetTolerance() const
    {
        return this->GetValue(AugmentedLagrangeConstraintAssembler::GetAlgorithmicParametersVariable())[1];
    }

    /// @internal
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
