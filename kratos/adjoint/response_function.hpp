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
#include "includes/kratos_export_api.h"
#include "includes/smart_pointers.h"

// --- STL Includes ---
#include <span>
#include <vector>


namespace Kratos {


/// @brief Class representing a differentiable response function.
/// @details The response function serves as an objective function,
///          or a component of another objective function. An important
///          restriction on the formulation of a reponse function, is that
///          it must be composable as a sum of elemental contributions.
///          @f[
///             j = \sum_e j_e
///          @f]
///          where
///          - @f$j@f$ is the response function,
///          - @f$e@f$ is an element index, and
///          - @f$j_e@f$ is the contribution of element @f$e@f$ to the response @f$j@f$.
///
///          This means that the response's derivatives must also be computable as a sum
///          of elemental contributions.
///          @f[
///             \frac{\partial j}{\partial \bullet} = \sum_e \frac{\partial j}{\partial \bullet}
///          @f]
/// @ingroup adjoints
class KRATOS_API(KRATOS_CORE) ResponseFunction : public IAdjoint {
public:
    KRATOS_CLASS_POINTER_DEFINITION(ResponseFunction);

    ResponseFunction() noexcept = default;

    /// @brief Construct a new instance and set its design variable types.
    /// @param[in] DesignVariableTypes List of design variables.
    /// @see @ref SetDesignVariableTypes
    ResponseFunction(std::span<const Variable<double>* const> DesignVariableTypes);

    /// @name Hooks
    /// @{

    /// @brief Function to be invoked before any calls to @p ComputeValue or @p ComputeDerivative in the current state.
    /// @param[in] rModelPart Computing @ref ModelPart.
    virtual void ComputeCache([[maybe_unused]] const ModelPart& rModelPart) {}

    /// @brief Function to be invoked after the last call to @p ComputeValue and @p ComputeDerivative in the current state.
    /// @param[in] rModelPart Computing @ref ModelPart.
    virtual void ClearCache() {}

    /// @}
    /// @name Evaluation
    /// @{

    /// @brief Compute the response value at the specified state.
    /// @param[in] rModelPart Computing @ref ModelPart.
    /// @param[in] iBuffer State variable index buffer indicating which step the response
    ///                    should be computed at, relative to the current one.
    [[nodiscard]] virtual double ComputeValue(
        const ModelPart& rModelPart,
        int iBuffer) const;

    /// @}
    /// @name Adjoint Interface
    /// @{

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

    /// @}

    /// @brief Specify which variables are to be treated as design variables.
    /// @param[in] Variables List of design variables.
    void SetDesignVariableTypes(std::span<const Variable<double>* const> Variables);

    /// @brief Fetch the stored set of design variable types.
    /// @param[out] rVariables Output list of design variable types.
    [[nodiscard]] std::span<const Variable<double>* const> GetDesignVariableTypes() const noexcept;

private:
    std::vector<const Variable<double>*> mDesignVariableTypes;
}; // class ResponseFunction


} // namespace Kratos
