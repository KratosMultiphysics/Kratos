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
    Element::Pointer pElement,
    const double Weight,
    const double ErrorThreshold)
        :   ResponseFunction(DesignVariableTypes),
            mName(rName),
            mpNode(pNode),
            mpElement(pElement),
            mWeight(Weight),
            mErrorThreshold(ErrorThreshold)
{}


Parameters SensorResponse::GetSensorParameters() const {
    KRATOS_TRY
        Parameters parameters = Parameters(R"(
        {
            "name"         : "",
            "location"     : [0.0, 0.0, 0.0],
            "weight"       : 0.0,
            "error_threshold": 0.0,
            "variable_data": {}
        })" );

        parameters["name"].SetString(this->GetName());
        parameters["location"].SetVector(this->GetNode()->Coordinates());
        parameters["weight"].SetDouble(this->GetWeight());
        parameters["error_threshold"].SetDouble(this->GetErrorThreshold());

        // Adding the data value container of the nodes
        const auto& r_node = *(this->GetNode());
        auto params_variable_data = parameters["variable_data"];
        for (const auto& r_var_value_pair : r_node.GetData()) {
            const std::string& var_name = std::get<0>(r_var_value_pair)->Name();
            if (KratosComponents<Variable<bool>>::Has(var_name)) {
                params_variable_data.AddBool(var_name, r_node.GetValue(KratosComponents<Variable<bool>>::Get(var_name)));
            } else if (KratosComponents<Variable<std::string>>::Has(var_name)) {
                params_variable_data.AddString(var_name, r_node.GetValue(KratosComponents<Variable<std::string>>::Get(var_name)));
            } else if (KratosComponents<Variable<int>>::Has(var_name)) {
                params_variable_data.AddInt(var_name, r_node.GetValue(KratosComponents<Variable<int>>::Get(var_name)));
            } else if (KratosComponents<Variable<double>>::Has(var_name)) {
                params_variable_data.AddDouble(var_name, r_node.GetValue(KratosComponents<Variable<double>>::Get(var_name)));
            } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(var_name)) {
                params_variable_data.AddVector(var_name, r_node.GetValue(KratosComponents<Variable<array_1d<double, 3>>>::Get(var_name)));
            } else if (KratosComponents<Variable<array_1d<double, 4>>>::Has(var_name)) {
                params_variable_data.AddVector(var_name, r_node.GetValue(KratosComponents<Variable<array_1d<double, 4>>>::Get(var_name)));
            } else if (KratosComponents<Variable<array_1d<double, 6>>>::Has(var_name)) {
                params_variable_data.AddVector(var_name, r_node.GetValue(KratosComponents<Variable<array_1d<double, 6>>>::Get(var_name)));
            } else if (KratosComponents<Variable<array_1d<double, 9>>>::Has(var_name)) {
                params_variable_data.AddVector(var_name, r_node.GetValue(KratosComponents<Variable<array_1d<double, 9>>>::Get(var_name)));
            } else if (KratosComponents<Variable<Vector>>::Has(var_name)) {
                params_variable_data.AddVector(var_name, r_node.GetValue(KratosComponents<Variable<Vector>>::Get(var_name)));
            } else if (KratosComponents<Variable<Matrix>>::Has(var_name)) {
                params_variable_data.AddMatrix(var_name, r_node.GetValue(KratosComponents<Variable<Matrix>>::Get(var_name)));
            } else {
                KRATOS_ERROR << "Unsupported variable type found under the variable name = " << var_name << ".\n";
            }
        }

        return parameters;
    KRATOS_CATCH("")
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
                << "sensor named \"" << this->GetName() << "\". Followings are available:";
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
    msg << "SensorResponse " << this->GetName();
    return msg.str();
}

void SensorResponse::PrintInfo(std::ostream& rOStream) const {
    rOStream << Info() << std::endl;
}

void SensorResponse::PrintData(std::ostream& rOStream) const {
    rOStream << "    Weight: " << this->GetWeight() << std::endl;
    rOStream << "    Error threshold: " << this->GetErrorThreshold() << std::endl;
    mpNode->PrintInfo(rOStream);
    mpNode->PrintData(rOStream);
}


} // namespace Kratos
