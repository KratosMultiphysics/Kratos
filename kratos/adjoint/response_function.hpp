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

// --- Kratos Core Includes ---
#include "adjoint/adjoint_interface.hpp"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/model_part.h"


namespace Kratos {


/// @brief Class representing a differentiable response function.
class KRATOS_API(KRATOS_CORE) ResponseFunction : public IAdjoint {
public:
    constexpr ResponseFunction() noexcept = default;

    /// @brief Construct a new instance and set its design variable types.
    /// @param[in] DesignVariableTypes List of design variables.
    /// @see @ref SetDesignVariableTypes
    ResponseFunction(std::span<const Variable<double>* const> DesignVariableTypes);

    /// @brief Compute the response value at the specified state.
    /// @param[in] rModelPart @ref ModelPart in the current state.
    /// @param[in] iBuffer State variable index buffer indicating which step the response
    ///                    should be computed at, relative to the current one.
    virtual double ComputeValue(
        const ModelPart& rModelPart,
        int iBuffer) const;

    /// @brief Specify which variables are to be treated as design variables.
    /// @param[in] Variables List of design variables.
    void SetDesignVariableTypes(std::span<const Variable<double>* const> Variables);

    /// @brief Fetch the stored set of design variable types.
    /// @param[out] rVariables Output list of design variable types.
    std::span<const Variable<double>* const> GetDesignVariableTypes() const noexcept;

    /// @brief Collect state variables from the provided element that influence the response function.
    /// @param[out] rVariables Output list of influencing state variables from @p rElement.
    /// @param[in] rElement @ref Element to collect variables from.
    /// @param[in] rProcessInfo Current @ref ProcessInfo of the computing @ref ModelPart.
    virtual void GetStateVariables(
        std::vector<IAdjoint::DynamicVariable>& rVariables,
        const Element& rElement,
        const ProcessInfo& rProcessInfo) const;

    /// @brief Collect state variables from the provided condition that influence the response function.
    /// @param[out] rVariables Output list of influencing state variables from @p rCondition.
    /// @param[in] rCondition @ref Condition to collect variables from.
    /// @param[in] rProcessInfo Current @ref ProcessInfo of the computing @ref ModelPart.
    virtual void GetStateVariables(
        std::vector<IAdjoint::DynamicVariable>& rVariables,
        const Condition& rCondition,
        const ProcessInfo& rProcessInfo) const;

    /// @brief Collect design variables from the provided element that influence the response function.
    /// @param[out] rVariables Output list of influencing design variables from @p rElement.
    /// @param[in] rElement @ref Element to collect variables from.
    /// @param[in] rProcessInfo Current @ref ProcessInfo of the computing @ref ModelPart.
    virtual void GetDesignVariables(
        std::vector<IAdjoint::DynamicVariable>& rVariables,
        const Element& rElement,
        const ProcessInfo& rProcessInfo) const;

    /// @brief Collect design variables from the provided condition that influence the response function.
    /// @param[out] rVariables Output list of influencing design variables from @p rCondition.
    /// @param[in] rCondition @ref Condition to collect variables from.
    /// @param[in] rProcessInfo Current @ref ProcessInfo of the computing @ref ModelPart.
    virtual void GetDesignVariables(
        std::vector<IAdjoint::DynamicVariable>& rVariables,
        const Condition& rCondition,
        const ProcessInfo& rProcessInfo) const;

    /// @brief Compute the response's derivative w.r.t the provided variables on the provided element.
    /// @param[out] rOutput Output vector containing the requested derivatives in the
    ///                     same order as @p Variables.
    /// @param[in] rElement Element providing the DoFs for computing the derivative.
    /// @param[in] Variables List of variables to compute the derivatives with respect to.
    /// @param[in] rProcessInfo Current @p ProcessInfo of the computing @ref ModelPart.
    /// @param[in] iBuffer State variable index buffer indicating which step the derivatives
    ///                    should be computed at, relative to the current one.
    virtual void ComputeDerivative(
        Vector& rOutput,
        const Element& rElement,
        std::span<const IAdjoint::DynamicVariable> Variables,
        const ProcessInfo& rProcessInfo,
        int iBuffer) const;

    /// @brief Compute the response's derivative w.r.t the provided variables on the provided condition.
    /// @param[out] rOutput Output vector containing the requested derivatives in the
    ///                     same order as @p Variables.
    /// @param[in] rCondition Condition providing the DoFs for computing the derivative.
    /// @param[in] Variables List of variables to compute the derivatives with respect to.
    /// @param[in] rProcessInfo Current @p ProcessInfo of the computing @ref ModelPart.
    /// @param[in] iBuffer State variable index buffer indicating which step the derivatives
    ///                    should be computed at, relative to the current one.
    virtual void ComputeDerivative(
        Vector& rOutput,
        const Condition& rCondition,
        std::span<const IAdjoint::DynamicVariable> Variables,
        const ProcessInfo& rProcessInfo,
        int iBuffer) const;

private:
    std::vector<const Variable<double>*> mDesignVariableTypes;
}; // class ResponseFunction


} // namespace Kratos
