//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:
//

#pragma once

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"
#include "containers/array_1d.h"
#include "containers/variable.h"
#include "includes/adjoint_interface.hpp"


namespace Kratos {


/// A base class for adjoint response functions.
class AdjointResponseFunction : public IAdjoint {
public:
    KRATOS_CLASS_POINTER_DEFINITION(AdjointResponseFunction);

    AdjointResponseFunction() noexcept = default;

    virtual ~AdjointResponseFunction() = default;

    ///@name Operations
    ///@{

    virtual void Initialize() {}

    virtual void InitializeSolutionStep() {}

    virtual void FinalizeSolutionStep() {}

    /// @brief Collect state variables from the provided element that influence the response function.
    /// @param[out] rVariables Output list of influencing state variables from @p rElement.
    /// @param[in] rElement @ref Element to collect variables from.
    /// @param[in] rProcessInfo Current @ref ProcessInfo of the computing @ref ModelPart.
    virtual void GetStateVariables(
        std::vector<IAdjoint::DynamicVariable>& rVariables,
        const Element& rElement,
        const ProcessInfo& rProcessInfo) const {
            KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
    }

    /// @brief Collect state variables from the provided condition that influence the response function.
    /// @param[out] rVariables Output list of influencing state variables from @p rCondition.
    /// @param[in] rCondition @ref Condition to collect variables from.
    /// @param[in] rProcessInfo Current @ref ProcessInfo of the computing @ref ModelPart.
    virtual void GetStateVariables(
        std::vector<IAdjoint::DynamicVariable>& rVariables,
        const Condition& rCondition,
        const ProcessInfo& rProcessInfo) const {
            KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
    }

    /// @brief Collect design variables from the provided element that influence the response function.
    /// @param[out] rVariables Output list of influencing design variables from @p rElement.
    /// @param[in] rElement @ref Element to collect variables from.
    /// @param[in] rProcessInfo Current @ref ProcessInfo of the computing @ref ModelPart.
    virtual void GetDesignVariables(
        std::vector<IAdjoint::DynamicVariable>& rVariables,
        const Element& rElement,
        const ProcessInfo& rProcessInfo) const {
            KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
    }

    /// @brief Collect design variables from the provided condition that influence the response function.
    /// @param[out] rVariables Output list of influencing design variables from @p rCondition.
    /// @param[in] rCondition @ref Condition to collect variables from.
    /// @param[in] rProcessInfo Current @ref ProcessInfo of the computing @ref ModelPart.
    virtual void GetDesignVariables(
        std::vector<IAdjoint::DynamicVariable>& rVariables,
        const Condition& rCondition,
        const ProcessInfo& rProcessInfo) const {
            KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
    }

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
        int iBuffer = 0) const {
            KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
    }

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
        int iBuffer = 0) const {
            KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
    }

    /// Calculate the local gradient w.r.t. primal solution.
    /**
     * @param[in]     rAdjointElement    the adjoint element.
     * @param[in]     rResidualGradient  the transposed gradient of the
     *                                   element's residual w.r.t. primal.
     * @param[out]    rResponseGradient  the gradient of the response function.
     * @param[in]     rProcessInfo       the current process info.
     */
    virtual void CalculateGradient(const Element& rAdjointElement,
                                   const Matrix& rResidualGradient,
                                   Vector& rResponseGradient,
                                   const ProcessInfo& rProcessInfo) {
        KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
    }

    /// Calculate the local gradient w.r.t. primal solution.
    /**
     * @param[in]     rAdjointCondition  the adjoint condition.
     * @param[in]     rResidualGradient  the transposed gradient of the
     *                                   condition's residual w.r.t. primal.
     * @param[out]    rResponseGradient  the gradient of the response function.
     * @param[in]     rProcessInfo       the current process info.
     */
    virtual void CalculateGradient(const Condition& rAdjointCondition,
                                   const Matrix& rResidualGradient,
                                   Vector& rResponseGradient,
                                   const ProcessInfo& rProcessInfo) {
        KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
    }

    /// Calculate the local gradient w.r.t. first derivatives of primal solution.
    /**
     * @param[in]     rAdjointElement    the adjoint element.
     * @param[in]     rResidualGradient  the transposed gradient of the
     *                                   element's residual w.r.t. first derivatives.
     * @param[out]    rResponseGradient  the gradient of the response function.
     * @param[in]     rProcessInfo       the current process info.
     */
    virtual void CalculateFirstDerivativesGradient(const Element& rAdjointElement,
                                                   const Matrix& rResidualGradient,
                                                   Vector& rResponseGradient,
                                                   const ProcessInfo& rProcessInfo) {
        KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
    }

    /// Calculate the local gradient w.r.t. first derivatives of primal solution.
    /**
     * @param[in]     rAdjointCondition  the adjoint condition.
     * @param[in]     rResidualGradient  the transposed gradient of the
     *                                   condition's residual w.r.t. first derivatives.
     * @param[out]    rResponseGradient  the gradient of the response function.
     * @param[in]     rProcessInfo       the current process info.
     */
    virtual void CalculateFirstDerivativesGradient(const Condition& rAdjointCondition,
                                                   const Matrix& rResidualGradient,
                                                   Vector& rResponseGradient,
                                                   const ProcessInfo& rProcessInfo) {
        KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
    }

    /// Calculate the local gradient w.r.t. second derivatives of primal solution.
    /**
     * @param[in]     rAdjointElement    the adjoint element.
     * @param[in]     rResidualGradient  the transposed gradient of the
     *                                   element's residual w.r.t. second derivatives.
     * @param[out]    rResponseGradient  the gradient of the response function.
     * @param[in]     rProcessInfo       the current process info.
     */
    virtual void CalculateSecondDerivativesGradient(const Element& rAdjointElement,
                                                    const Matrix& rResidualGradient,
                                                    Vector& rResponseGradient,
                                                    const ProcessInfo& rProcessInfo) {
        KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
    }

    /// Calculate the local gradient w.r.t. second derivatives of primal solution.
    /**
     * @param[in]     rAdjointCondition  the adjoint condition.
     * @param[in]     rResidualGradient  the transposed gradient of the
     *                                   condition's residual w.r.t. second derivatives.
     * @param[out]    rResponseGradient  the gradient of the response function.
     * @param[in]     rProcessInfo       the current process info.
     */
    virtual void CalculateSecondDerivativesGradient(const Condition& rAdjointCondition,
                                                    const Matrix& rResidualGradient,
                                                    Vector& rResponseGradient,
                                                    const ProcessInfo& rProcessInfo) {
        KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
    }

    /// Calculate the partial sensitivity w.r.t. design variable.
    /**
     * @param[in]     rAdjointElement       the adjoint element.
     * @param[in]     rVariable             the design variable.
     * @param[in]     rSensitivityMatrix    the transposed gradient of the
     *                                      element's residual w.r.t. design variable.
     * @param[out]    rSensitivityGradient  the gradient of the response function.
     * @param[in]     rProcessInfo          the current process info.
     */
    virtual void CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo) {
        KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
    }

    /// Calculate the partial sensitivity w.r.t. design variable.
    /**
     * @param[in]     rAdjointCondition     the adjoint condition.
     * @param[in]     rVariable             the design variable.
     * @param[in]     rSensitivityMatrix    the transposed gradient of the
     *                                      condition's residual w.r.t. design variable.
     * @param[out]    rSensitivityGradient  the gradient of the response function.
     * @param[in]     rProcessInfo          the current process info.
     */
    virtual void CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo) {
        KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
    }

    /// Calculate the partial sensitivity w.r.t. design variable.
    /**
     * @param[in]     rAdjointElement       the adjoint element.
     * @param[in]     rVariable             the design variable.
     * @param[in]     rSensitivityMatrix    the transposed gradient of the
     *                                      element's residual w.r.t. design variable.
     * @param[out]    rSensitivityGradient  the gradient of the response function.
     * @param[in]     rProcessInfo          the current process info.
     */
    virtual void CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo) {
        KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
    }

    /// Calculate the partial sensitivity w.r.t. design variable.
    /**
     * @param[in]     rAdjointCondition     the adjoint condition.
     * @param[in]     rVariable             the design variable.
     * @param[in]     rSensitivityMatrix    the transposed gradient of the
     *                                      condition's residual w.r.t. design variable.
     * @param[out]    rSensitivityGradient  the gradient of the response function.
     * @param[in]     rProcessInfo          the current process info.
     */
    virtual void CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo) {
        KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
    }

    /// Calculate the scalar valued response function.
    virtual double CalculateValue(ModelPart& rModelPart) = 0;
}; // class AdjointResponseFunction

} // namespace Kratos
