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
#include "adjoint/response_function.hpp" // ResponseFunction
#include "includes/code_location.h" // KRATOS_CODE_LOCATION


namespace Kratos {


ResponseFunction::ResponseFunction(std::span<const Variable<double>* const> DesignVariableTypes) {
    this->SetDesignVariableTypes(DesignVariableTypes);
}


double ResponseFunction::ComputeValue(
    const Element&,
    const ProcessInfo&,
    int) const {
        KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
}


double ResponseFunction::ComputeValue(
    const Condition&,
    const ProcessInfo&,
    int) const {
        KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
}


void ResponseFunction::SetDesignVariableTypes(std::span<const Variable<double>* const> Variables) {
    KRATOS_TRY
        mDesignVariableTypes.clear();
        mDesignVariableTypes.insert(
            mDesignVariableTypes.end(),
            Variables.begin(),
            Variables.end());
        std::sort(
            mDesignVariableTypes.begin(),
            mDesignVariableTypes.end(),
            [] (const auto p_left, const auto p_right) -> bool {
                return p_left->Key() < p_right->Key();
            });
    KRATOS_CATCH("")
}


std::span<const Variable<double>* const> ResponseFunction::GetDesignVariableTypes() const noexcept {
    return mDesignVariableTypes;
}


void ResponseFunction::GetStateVariables(
    std::vector<IAdjoint::DynamicVariable>&,
    const Element&,
    const ProcessInfo&) const {
        KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
}


void ResponseFunction::GetStateVariables(
    std::vector<IAdjoint::DynamicVariable>&,
    const Condition&,
    const ProcessInfo&) const {
        KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
}


/// @details Remove entries from @p rVariables whose keys are not
///          present in @p DesignVariables. @p DesignVariables are
///          assumed to be ordered with respect to their keys.
void VariableSetIntersection(
    std::vector<IAdjoint::DynamicVariable>& rVariables,
    std::span<const Variable<double>* const> DesignVariables) {
        KRATOS_TRY
            rVariables.erase(
                std::remove_if(
                    rVariables.begin(),
                    rVariables.end(),
                    [&] (const IAdjoint::DynamicVariable& r_variable) -> bool {
                        const auto it_design_variable = std::lower_bound(
                            DesignVariables.begin(),
                            DesignVariables.end(),
                            r_variable,
                            [] (const Variable<double>* p_left, const IAdjoint::DynamicVariable& r_right) -> bool {
                                return p_left->Key() < r_right.Key();
                            });
                        return it_design_variable == DesignVariables.end()
                            || (*it_design_variable)->Key() != r_variable.Key();
                    }),
                rVariables.end());
        KRATOS_CATCH("")
}


void ResponseFunction::GetDesignVariables(
    std::vector<IAdjoint::DynamicVariable>& rVariables,
    const Element& rElement,
    const ProcessInfo& rProcessInfo) const {
        KRATOS_TRY
            rElement.GetInfluencingVariables(
                rVariables,
                rProcessInfo);
            VariableSetIntersection(
                rVariables,
                mDesignVariableTypes);
        KRATOS_CATCH("")
}


void ResponseFunction::GetDesignVariables(
    std::vector<IAdjoint::DynamicVariable>& rVariables,
    const Condition& rCondition,
    const ProcessInfo& rProcessInfo) const {
        KRATOS_TRY
            rCondition.GetInfluencingVariables(
                rVariables,
                rProcessInfo);
            VariableSetIntersection(
                rVariables,
                mDesignVariableTypes);
        KRATOS_CATCH("")
}


void ResponseFunction::ComputeDerivative(
    Vector& rOutput,
    const Element& rElement,
    std::span<const IAdjoint::DynamicVariable> Variables,
    const ProcessInfo& rProcessInfo,
    int iBuffer) const {
        KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
}


void ResponseFunction::ComputeDerivative(
    Vector& rOutput,
    const Condition& rCondition,
    std::span<const IAdjoint::DynamicVariable> Variables,
    const ProcessInfo& rProcessInfo,
    int iBuffer) const {
        KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
}


} // namespace Kratos
