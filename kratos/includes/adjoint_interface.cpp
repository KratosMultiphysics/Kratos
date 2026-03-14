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
#include "includes/adjoint_interface.hpp" // IAdjoint, IAdjointElement
#include "includes/exception.h" // KRATOS_ERROR
#include "includes/code_location.h" // KRATOS_CODE_LOCATION

// --- STL Includes ---
#include <algorithm> // std::unique, std::sort


namespace Kratos {


void IAdjoint::GetStateVariables(std::vector<VARIABLE>&) const {
    KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
}


void IAdjoint::GetInfluencingVariables(std::vector<IAdjoint::VARIABLE>&) const {
    KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
}


void IAdjointElement::GetDofs(std::vector<const Dof<IAdjoint::Scalar>*>&) const {
    KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
}


void IAdjointElement::GetInfluencingVariables(std::vector<IAdjoint::VARIABLE>& rOutput) const {
    KRATOS_TRY
        std::vector<IAdjoint::VARIABLE> buffer;
        #define KRATOS_IADJOINTELEMENT_COLLECT(term)                    \
            this->GetInfluencingVariables<term>(buffer);                \
            rOutput.insert(rOutput.end(), buffer.begin(), buffer.end())
        KRATOS_IADJOINTELEMENT_COLLECT(IAdjoint::ResidualTerm::Mass);
        KRATOS_IADJOINTELEMENT_COLLECT(IAdjoint::ResidualTerm::Damping);
        KRATOS_IADJOINTELEMENT_COLLECT(IAdjoint::ResidualTerm::Stiffness);
        KRATOS_IADJOINTELEMENT_COLLECT(IAdjoint::ResidualTerm::Load);
        #undef KRATOS_IADJOINTELEMENT_COLLECT

        std::sort(
            rOutput.begin(),
            rOutput.end(),
            [] (const IAdjoint::VARIABLE& r_left, const IAdjoint::VARIABLE& r_right) -> bool {
                return r_left.Key() < r_right.Key();
            });
        rOutput.erase(
            std::unique(
                rOutput.begin(),
                rOutput.end(),
                [] (const IAdjoint::VARIABLE& r_left, const IAdjoint::VARIABLE& r_right) -> bool {
                    return r_left.Key() == r_right.Key();
                }),
            rOutput.end());
    KRATOS_CATCH("")
}


void IAdjointElement::GetMassInfluencingVariables(std::vector<IAdjoint::VARIABLE>&) const {
    KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
}


void IAdjointElement::GetDampingInfluencingVariables(std::vector<IAdjoint::VARIABLE>&) const {
    KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
}


void IAdjointElement::GetStiffnessInfluencingVariables(std::vector<IAdjoint::VARIABLE>&) const {
    KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
}


void IAdjointElement::GetLoadInfluencingVariables(std::vector<IAdjoint::VARIABLE>&) const {
    KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
}


void IAdjointElement::ComputeStiffnessDerivative(
    std::span<const IAdjoint::VARIABLE> Variables,
    Matrix& rOutput) const {
        KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
}


void IAdjointElement::ComputeDampingDerivative(
    std::span<const IAdjoint::VARIABLE> Variables,
    Matrix& rOutput) const {
        KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
}


void IAdjointElement::ComputeMassDerivative(
    std::span<const IAdjoint::VARIABLE> Variables,
    Matrix& rOutput) const {
        KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
}


void IAdjointElement::ComputeLoadDerivative(
    std::span<const IAdjoint::VARIABLE> Variables,
    Matrix& rOutput) const {
        KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
}


} // namespace Kratos
