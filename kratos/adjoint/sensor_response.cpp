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
#include "adjoint/sensor_response.hpp"


namespace Kratos {


SensorResponse::SensorResponse(
    std::span<const Variable<double>* const> DesignVariableTypes,
    const std::string& rName,
    Node::Pointer pNode,
    std::optional<intrusive_ptr<const Element>> pMaybeElement)
        :   ResponseFunction(DesignVariableTypes),
            mName(rName),
            mpNode(pNode),
            mpMaybeElement(pMaybeElement)
{}


void SensorResponse::GetStateVariables(
    std::vector<IAdjoint::DynamicVariable>& rVariables,
    const Condition&,
    const ProcessInfo&) const {
        rVariables.clear();
}


void SensorResponse::GetDesignVariables(
    std::vector<IAdjoint::DynamicVariable>& rVariables,
    const Condition&,
    const ProcessInfo&) const {
        rVariables.clear();
}


void SensorResponse::ComputeDerivative(
    Vector& rOutput,
    const Condition&,
    std::span<const IAdjoint::DynamicVariable> Variables,
    const ProcessInfo&,
    int) const {
        for (const auto& r_variable : Variables) {
            // Check whether the requested variable is a design variable.
            const auto& r_design_variable_types = this->GetDesignVariableTypes();
            const auto it_design_variable_type = std::find_if(
                r_design_variable_types.begin(),
                r_design_variable_types.end(),
                [&r_variable] (const Variable<double>* p_variable_type) -> bool {
                    return p_variable_type->Key() == r_variable;
                });
            const bool is_design_variable = it_design_variable_type != r_design_variable_types.end();

            // Error if the requested variable is neither a state nor a design variable.
            KRATOS_ERROR_IF_NOT(is_design_variable) << "unsupported variable " << r_variable.Name();
        } // for r_variable in Variables
        rOutput = ZeroVector(Variables.size());
}


void SensorResponse::AddTensorAdaptor(
    const std::string& rTensorAdaptorName,
    TensorAdaptor<double>::Pointer pTensorAdaptor) {
    KRATOS_TRY
        const auto p_itr = mTensorAdaptorsMap.find(rTensorAdaptorName);
        KRATOS_ERROR_IF_NOT(p_itr == mTensorAdaptorsMap.end())
            << "A tensor adaptor named \"" << rTensorAdaptorName << " already exists.";
        mTensorAdaptorsMap[rTensorAdaptorName] = pTensorAdaptor;
    KRATOS_CATCH("");
}

TensorAdaptor<double>::Pointer SensorResponse::GetTensorAdaptor(const std::string& rTensorAdaptorName) const {
    KRATOS_TRY
        const auto p_itr = mTensorAdaptorsMap.find(rTensorAdaptorName);

        if (p_itr == mTensorAdaptorsMap.end()) {
            std::stringstream msg;
            msg << "A tensor adaptor named \"" << rTensorAdaptorName << "\" not found in "
                << "sensor named \"" << this->Name() << "\". Followings are available:";
            for (const auto& r_pair : mTensorAdaptorsMap) {
                msg << std::endl << "   "  << r_pair.first;
            }
            KRATOS_ERROR << msg.str();
        }

        return p_itr->second;
    KRATOS_CATCH("");
}


void SensorResponse::ClearTensorAdaptors() {
    mTensorAdaptorsMap.clear();
}


std::string SensorResponse::Info() const {
    std::stringstream msg;
    msg << "SensorResponse " << this->Name();
    return msg.str();
}

void SensorResponse::PrintInfo(std::ostream& rOStream) const {
    rOStream << this->Info() << std::endl;
}

void SensorResponse::PrintData(std::ostream& rOStream) const {
    mpNode->PrintInfo(rOStream);
    mpNode->PrintData(rOStream);
}


} // namespace Kratos
