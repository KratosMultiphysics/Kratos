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
///          restriction on the formulation of a reponse is that
///          it must be expressable as
///          @f[
///             j = f \left( \sum_e j_e \right)
///          @f]
///          where
///          - @f$j@f$ is the response function,
///          - @f$f@f$ is the outer response (see @ref ComputeValue and @ref ComputeDerivative),
///          - @f$e@f$ is an element index, and
///          - @f$j_e@f$ is the inner response from element @f$e@f$ (see @ref ComputeInnerValue and @ref ComputeInnerDerivative).
///
///          This means that the response's derivatives take the following form:
///          @f[
///             \frac{\partial j}{\partial \bullet} = \frac{\partial f}{\partial \bullet} \sum_e \frac{\partial j_e}{\partial \bullet}
///          @f]
///
///          Pseudocode computing the response value:
///          @code{.py}
///             model_part: ModelPart
///             response: ResponseFunction
///
///             response.ComputeCache(model_part)
///             inner_response: float = 0.0
///             for element in elements:
///                 inner_response += response.ComputeInnerValue(
///                     element,
///                     model_part.ProcessInfo,
///                     0)
///
///             response_value = response.ComputeValue(inner_response)
///          @endcode
///
///          Pseudocode computing the response derivative (e.g.: w.r.t state variables):
///          @code{.py}
///             model_part: ModelPart
///             response: ResponseFunction
///
///             response.ComputeCache(model_part)
///             inner_derivative: numpy.ndarray = numpy.zeros(system_size)
///             for element in elements:
///                 state_variables = response.GetStateVariables(
///                     element,
///                     model_part.ProcessInfo)
///                 global_indices = element.EquationIdVector(model_part.ProcessInfo)
///                 element_inner_derivative = response.ComputeInnerDerivative(
///                     element,
///                     state_variables,
///                     model_part.ProcessInfo,
///                     0)
///                 # ... assemble element_inner_response into inner_derivative
///
///             response_derivative = response.ComputeDerivative(inner_derivative)
///          @endcode
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

    /// @brief Function to be invoked before any calls to @p ComputeValue, @p ComputeInnerDerivative or @p ComputeDerivative in the current state.
    /// @param[in] rModelPart Computing @ref ModelPart.
    virtual void ComputeCache([[maybe_unused]] const ModelPart& rModelPart) {}

    /// @brief Function to be invoked after the last call to @p ComputeValue, @p ComputeInnerDerivative and @p ComputeDerivative in the current state.
    /// @param[in] rModelPart Computing @ref ModelPart.
    virtual void ClearCache() {}

    /// @}
    /// @name Evaluation
    /// @{

    /// @brief Compute @f$j_e@f$, a component of the inner response value at the specified state.
    /// @param[in] rElement Element to compute contributions for.
    /// @param[in] rProcessInfo Current @ref ProcessInfo of the computing @ref ModelPart.
    /// @param[in] iBuffer State variable index buffer indicating which step the response
    ///                    should be computed at, relative to the current one.
    [[nodiscard]] virtual double ComputeInnerValue(
        const Element& rElement,
        const ProcessInfo& rProcessInfo,
        int iBuffer) const;

    /// @brief Compute @f$j_e@f$, a component of the inner response value at the specified state.
    /// @param[in] rCondition Condition to compute contributions for.
    /// @param[in] rProcessInfo Current @ref ProcessInfo of the computing @ref ModelPart.
    /// @param[in] iBuffer State variable index buffer indicating which step the response
    ///                    should be computed at, relative to the current one.
    [[nodiscard]] virtual double ComputeInnerValue(
        const Condition& rCondition,
        const ProcessInfo& rProcessInfo,
        int iBuffer) const;

    /// @brief Compute the response @f$j@f$, given the inner response @f$\sum_e j_e@f$.
    /// @param[in] InnerSum Summed inner response value @f$\sum_e j_e@f$.
    /// @param[in] rProcessInfo Current @ref ProcessInfo of the computing @ref ModelPart.
    /// @param[in] iBuffer State variable index buffer indicating which step the response
    ///                    should be computed at, relative to the current one.
    [[nodiscard]] virtual double ComputeValue(
        double InnerSum,
        const ProcessInfo& rProcessInfo,
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

    /// @brief Compute the derivative of @f$j_e@f$ w.r.t the specified variables on the provided element.
    /// @param[out] rOutput Output vector containing the requested derivatives in the
    ///                     same order as @p Variables.
    /// @param[in] rElement Element providing the DoFs for computing the derivative.
    /// @param[in] Variables List of variables to compute the derivatives with respect to.
    /// @param[in] rProcessInfo Current @p ProcessInfo of the computing @ref ModelPart.
    /// @param[in] iBuffer State variable index buffer indicating which step the derivatives
    ///                    should be computed at, relative to the current one.
    virtual void ComputeInnerDerivative(
        Vector& rOutput,
        const Element& rElement,
        std::span<const IAdjoint::DynamicVariable> Variables,
        const ProcessInfo& rProcessInfo,
        int iBuffer) const;

    /// @brief Compute the derivative of @f$j_e@f$ w.r.t the specified variables on the provided condition.
    /// @param[out] rOutput Output vector containing the requested derivatives in the
    ///                     same order as @p Variables.
    /// @param[in] rCondition Condition providing the DoFs for computing the derivative.
    /// @param[in] Variables List of variables to compute the derivatives with respect to.
    /// @param[in] rProcessInfo Current @p ProcessInfo of the computing @ref ModelPart.
    /// @param[in] iBuffer State variable index buffer indicating which step the derivatives
    ///                    should be computed at, relative to the current one.
    virtual void ComputeInnerDerivative(
        Vector& rOutput,
        const Condition& rCondition,
        std::span<const IAdjoint::DynamicVariable> Variables,
        const ProcessInfo& rProcessInfo,
        int iBuffer) const;

    /// @brief Compute @f$\frac{\partial j}{\partial \bullet}@f$, given @f$\sum_e \frac{\partial j_e}{\partial \bullet}@f$.
    /// @param[inout] rDerivative Passed as the inner sum computed in @ref ComputeInnerDerivative, and modified by
    ///                           this function to become @f$\frac{\partial j}{\partial \bullet}@f$.
    /// @param[in] Variables List of variables to compute the derivatives with respect to.
    /// @param[in] rProcessInfo Current @p ProcessInfo of the computing @ref ModelPart.
    /// @param[in] iBuffer State variable index buffer indicating which step the derivatives
    ///                    should be computed at, relative to the current one.
    virtual void ComputeDerivative(
        Vector& rDerivative,
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
