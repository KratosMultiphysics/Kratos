//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         SystemIdentificationApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <type_traits>

// External includes

// Project includes

// Application includes
#include "includes/kratos_components.h"

// Include base h
#include "sensor.h"

namespace Kratos
{

Sensor::Sensor(
    const std::string& rName,
    Node::Pointer pNode,
    const double Weight)
    : mName(rName),
      mpNode(pNode),
      mWeight(Weight),
      mSensorValue(0.0)
{
}

Parameters Sensor::GetSensorParameters() const
{
    Parameters parameters = Parameters(R"(
    {
        "name"         : "",
        "value"        : 0.0,
        "location"     : [0.0, 0.0, 0.0],
        "weight"       : 0.0,
        "variable_data": {}
    })" );

    parameters["name"].SetString(this->GetName());
    parameters["value"].SetDouble(this->GetSensorValue());
    parameters["location"].SetVector(this->GetNode()->Coordinates());
    parameters["weight"].SetDouble(this->GetWeight());

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
}

std::string Sensor::GetName() const
{
    return mName;
}

Node::Pointer Sensor::GetNode() const
{
    return mpNode;
}

double Sensor::GetWeight() const
{
    return mWeight;
}

double Sensor::GetSensorValue() const
{
    return mSensorValue;
}

void Sensor::SetSensorValue(const double Value)
{
    mSensorValue = Value;
}

void Sensor::AddTensorAdaptor(
    const std::string& rTensorAdaptorName,
    TensorAdaptor<double>::Pointer pTensorAdaptor)
{
    KRATOS_TRY

    const auto p_itr = mTensorAdaptorsMap.find(rTensorAdaptorName);
    KRATOS_ERROR_IF_NOT(p_itr == mTensorAdaptorsMap.end())
        << "A tensor adaptor named \"" << rTensorAdaptorName << " already exists.";
    mTensorAdaptorsMap[rTensorAdaptorName] = pTensorAdaptor;

    KRATOS_CATCH("");
}

TensorAdaptor<double>::Pointer Sensor::GetTensorAdaptor(const std::string& rTensorAdaptorName) const
{
    KRATOS_TRY

    const auto p_itr = mTensorAdaptorsMap.find(rTensorAdaptorName);

    if (p_itr == mTensorAdaptorsMap.end()) {
        std::stringstream msg;
        msg << "A tensor adaptor named \"" << rTensorAdaptorName << "\" not found in "
            << "sensor named \"" << this->GetName() << "\". Followings are available:";
        for (const auto& r_pair : this->GetTensorAdaptorsMap()) {
            msg << std::endl << "   "  << r_pair.first;
        }
        KRATOS_ERROR << msg.str();
    }

    return p_itr->second;

    KRATOS_CATCH("");
}

std::unordered_map<std::string, TensorAdaptor<double>::Pointer> Sensor::GetTensorAdaptorsMap() const
{
    return mTensorAdaptorsMap;
}

void Sensor::ClearTensorAdaptors()
{
    mTensorAdaptorsMap.clear();
}

std::string Sensor::Info() const
{
    std::stringstream msg;
    msg << "Sensor " << this->GetName();
    return msg.str();
}

void Sensor::PrintInfo(std::ostream& rOStream) const
{
    rOStream << Info() << std::endl;
}

void Sensor::PrintData(std::ostream& rOStream) const
{
    rOStream << "    Value: " << this->GetSensorValue() << std::endl;
    rOStream << "    Weight: " << this->GetWeight() << std::endl;
    mpNode->PrintInfo(rOStream);
    mpNode->PrintData(rOStream);
}

} /* namespace Kratos.*/