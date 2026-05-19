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

// --- Kratos Includes ---
#include "includes/adjoint_interface.hpp"
#include "includes/element.h"
#include "includes/condition.h"

// --- STL Includes ---
#include <span> // std::span



namespace Kratos {


template <IAdjoint::ResidualTerm Term, class TEntity>
requires (std::is_same_v<TEntity,Element> || std::is_same_v<TEntity,Condition>)
class AdjointFiniteDifferenceUtility {
public:
    /// @brief Class uniquely defining what to perturb in an element.
    struct Perturbation {
        /// @brief Variable to perturb.
        /// @details @ref IAdjoint::DynamicVariable optionally offers to augment the variable
        ///          with an integer whose interpretation depends on the @ref mContext "specified context".
        ///          For example, in the context of an element and a @p YOUNG_MODULUS variable, it can refer
        ///          to the Young's modulus of a specific layer within a multi-layered composite material.
        IAdjoint::DynamicVariable mVariable;

        /// @brief Context specifying which container the @ref mVariable "variable" refers to.
        Globals::DataLocation mContext;

        /// @brief Magnitude of the perturbation to apply.
        double mMagnitude;
    }; // struct Perturbation

    using Value = std::conditional_t<
        Term == IAdjoint::ResidualTerm::Load,
        Vector,
        Matrix>;

    void FiniteDifferenceDerivative(
        const TEntity& rEntity,
        const Vector& rValue,
        std::span<const Perturbation> Perturbations,
        Matrix& rOutput,
        std::size_t iBuffer,
        const ProcessInfo& rProcessInfo) const
    requires (Term != IAdjoint::ResidualTerm::Load);

    void FiniteDifferenceDerivative(
        const TEntity& rEntity,
        std::span<const Perturbation> Perturbations,
        Matrix& rOutput,
        std::size_t iBuffer,
        const ProcessInfo& rProcessInfo) const
    requires (Term == IAdjoint::ResidualTerm::Load);

private:
    void FiniteDifferenceDerivative(
        const TEntity& rEntity,
        const Vector* pValue,
        std::span<const Perturbation> Perturbations,
        Matrix& rOutput,
        std::size_t iBuffer,
        const ProcessInfo& rProcessInfo) const;

    static void EvaluateTerm(
        TEntity& rEntity,
        const ProcessInfo& rProcessInfo,
        Value& rOutput);

    template <Globals::DataLocation TContext>
    void Perturb(
        TEntity& rEntity,
        const IAdjoint::DynamicVariable& rVariable,
        double PerturbationMagnitude,
        std::size_t iBuffer) const;

    void ComputeFiniteDifferences(
        Value& rOutput,
        const Value& rReferenceState,
        const Value& rPerturbedState,
        double PerturbationMagnitude,
        std::size_t OutputBegin) const;

    static void ResizeOutput(
        Matrix& rOutput,
        std::span<const std::size_t> SingleVariateDerivativeShape,
        std::size_t VariableCount);
}; // class AdjointFiniteDifferenceUtility


} // namespace Kratos
