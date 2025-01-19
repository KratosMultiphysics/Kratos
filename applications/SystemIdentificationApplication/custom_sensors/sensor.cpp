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
        "name"      : "",
        "value"     : 0.0,
        "location"  : [0.0, 0.0, 0.0],
        "weight"    : 0.0
    })" );

    parameters["name"].SetString(this->GetName());
    parameters["value"].SetDouble(this->GetSensorValue());
    parameters["location"].SetVector(this->GetNode()->Coordinates());
    parameters["weight"].SetDouble(this->GetWeight());

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

void Sensor::AddContainerExpression(
    const std::string& rExpressionName,
    ContainerExpressionType pContainerExpression)
{
    KRATOS_TRY

    const auto p_itr = mContainerExpressions.find(rExpressionName);
    KRATOS_ERROR_IF_NOT(p_itr == mContainerExpressions.end())
        << "A container expression named \"" << rExpressionName << " already exists.";
    mContainerExpressions[rExpressionName] = pContainerExpression;

    KRATOS_CATCH("");
}

Sensor::ContainerExpressionType Sensor::GetContainerExpression(const std::string& rExpressionName) const
{
    KRATOS_TRY

    const auto p_itr = mContainerExpressions.find(rExpressionName);

    KRATOS_ERROR_IF(p_itr == mContainerExpressions.end())
        << "A container expression named \"" << rExpressionName << " not found.";

    return p_itr->second;

    KRATOS_CATCH("");
}

std::unordered_map<std::string, Sensor::ContainerExpressionType> Sensor::GetContainerExpressionsMap() const
{
    return mContainerExpressions;
}

void Sensor::ClearContainerExpressions()
{
    mContainerExpressions.clear();
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