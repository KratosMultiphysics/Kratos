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

// --- Kratos Core Includes ---
#include "includes/condition.h"


namespace Kratos {


void Condition::GetMassInfluencingVariables(
    std::vector<IAdjoint::DynamicVariable>& rOutput,
    const ProcessInfo&) const {
        rOutput.clear();
}


void Condition::GetDampingInfluencingVariables(
    std::vector<IAdjoint::DynamicVariable>& rOutput,
    const ProcessInfo&) const {
        rOutput.clear();
}


void Condition::GetStiffnessInfluencingVariables(
    std::vector<IAdjoint::DynamicVariable>& rOutput,
    const ProcessInfo&) const {
        rOutput.clear();
}


void Condition::ComputeMassDerivative(
    Matrix& rOutput,
    std::span<const IAdjoint::DynamicVariable> Variables,
    const ProcessInfo&,
    int iBuffer) const {
        rOutput = ZeroMatrix(Variables.size(), Variables.size());
}


void Condition::ComputeDampingDerivative(
    Matrix& rOutput,
    std::span<const IAdjoint::DynamicVariable> Variables,
    const ProcessInfo&,
    int iBuffer) const {
        rOutput = ZeroMatrix(Variables.size(), Variables.size());
}


void Condition::ComputeStiffnessDerivative(
    Matrix& rOutput,
    std::span<const IAdjoint::DynamicVariable> Variables,
    const ProcessInfo&,
    int iBuffer) const {
        rOutput = ZeroMatrix(Variables.size(), Variables.size());
}


} // namespace Kratos
