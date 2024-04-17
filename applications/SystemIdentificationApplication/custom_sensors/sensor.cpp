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
    const Point& rLocation,
    const double Weight)
    : mName(rName),
      mLocation(rLocation),
      mWeight(Weight),
      mSensorValue(0.0)
{
}

std::string Sensor::GetName() const
{
    return mName;
}

Point Sensor::GetLocation() const
{
    return mLocation;
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

template<class TContainerType>
void Sensor::AddContainerExpression(
    const std::string& rExpressionName,
    typename ContainerExpression<TContainerType>::Pointer pContainerExpression)
{
    KRATOS_TRY

    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        const auto p_itr = mNodalExpressions.find(rExpressionName);
        KRATOS_ERROR_IF_NOT(p_itr == mNodalExpressions.end())
            << "A nodal expression named \"" << rExpressionName << " already exists.";
        mNodalExpressions[rExpressionName] = pContainerExpression;
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
        const auto p_itr = mConditionExpressions.find(rExpressionName);
        KRATOS_ERROR_IF_NOT(p_itr == mConditionExpressions.end())
            << "A condition expression named \"" << rExpressionName << " already exists.";
        mConditionExpressions[rExpressionName] = pContainerExpression;
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {
        const auto p_itr = mElementExpressions.find(rExpressionName);
        KRATOS_ERROR_IF_NOT(p_itr == mElementExpressions.end())
            << "A element expression named \"" << rExpressionName << " already exists.";
        mElementExpressions[rExpressionName] = pContainerExpression;
    } else {
        static_assert(std::is_same_v<TContainerType, TContainerType>, "Unsupported container type.");
    }

    KRATOS_CATCH("");
}

void Sensor::AddNodalExpression(
    const std::string& rExpressionName,
    ContainerExpression<ModelPart::NodesContainerType>::Pointer pNodalExpression)
{
    AddContainerExpression<ModelPart::NodesContainerType>(rExpressionName, pNodalExpression);
}

ContainerExpression<ModelPart::NodesContainerType>::Pointer Sensor::GetNodalExpression(const std::string& rExpressionName) const
{
    KRATOS_TRY

    const auto p_itr = mNodalExpressions.find(rExpressionName);

    KRATOS_ERROR_IF(p_itr == mNodalExpressions.end())
        << "A nodal expression named \"" << rExpressionName << " not found.";

    return p_itr->second;

    KRATOS_CATCH("");
}

std::unordered_map<std::string, ContainerExpression<ModelPart::NodesContainerType>::Pointer> Sensor::GetNodalExpressionsMap() const
{
    return mNodalExpressions;
}

void Sensor::AddConditionExpression(
    const std::string& rExpressionName,
    ContainerExpression<ModelPart::ConditionsContainerType>::Pointer pConditionExpression)
{
    AddContainerExpression<ModelPart::ConditionsContainerType>(rExpressionName, pConditionExpression);
}

ContainerExpression<ModelPart::ConditionsContainerType>::Pointer Sensor::GetConditionExpression(const std::string& rExpressionName) const
{
    KRATOS_TRY

    const auto p_itr = mConditionExpressions.find(rExpressionName);

    KRATOS_ERROR_IF(p_itr == mConditionExpressions.end())
        << "A condition expression named \"" << rExpressionName << " not found.";

    return p_itr->second;

    KRATOS_CATCH("");
}

std::unordered_map<std::string, ContainerExpression<ModelPart::ConditionsContainerType>::Pointer> Sensor::GetConditionExpressionsMap() const
{
    return mConditionExpressions;
}

void Sensor::AddElementExpression(
    const std::string& rExpressionName,
    ContainerExpression<ModelPart::ElementsContainerType>::Pointer pElementExpression)
{
    AddContainerExpression<ModelPart::ElementsContainerType>(rExpressionName, pElementExpression);
}

ContainerExpression<ModelPart::ElementsContainerType>::Pointer Sensor::GetElementExpression(const std::string& rExpressionName) const
{
    KRATOS_TRY

    const auto p_itr = mElementExpressions.find(rExpressionName);

    KRATOS_ERROR_IF(p_itr == mElementExpressions.end())
        << "A condition expression named \"" << rExpressionName << " not found.";

    return p_itr->second;

    KRATOS_CATCH("");
}

std::unordered_map<std::string, ContainerExpression<ModelPart::ElementsContainerType>::Pointer> Sensor::GetElementExpressionsMap() const
{
    return mElementExpressions;
}

std::vector<std::string> Sensor::GetDataVariableNames() const
{
    std::vector<std::string> result;
    for (auto i_ptr = this->begin(); i_ptr != this->end(); ++i_ptr) {
        result.push_back(i_ptr->first->Name());
    }
    return result;
}

void Sensor::ClearNodalExpressions()
{
    mNodalExpressions.clear();
}

void Sensor::ClearConditionExpressions()
{
    mConditionExpressions.clear();
}

void Sensor::ClearElementExpressions()
{
    mElementExpressions.clear();
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
    PrintInfo(rOStream);
    rOStream << "    Location: " << this->GetLocation() << std::endl;
    rOStream << "    Value: " << this->GetSensorValue() << std::endl;
    rOStream << "    Weight: " << this->GetWeight() << std::endl;
    DataValueContainer::PrintData(rOStream);
}

} /* namespace Kratos.*/