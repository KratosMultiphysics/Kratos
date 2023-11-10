//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: DigitalTwinApplication/license.txt
//
//  Main authors:    Suneth Wranakulasuriya
//

// System includes

// External includes

// Project includes

// Application includes

// Include base h
#include "sensor_specification.h"

namespace Kratos
{

SensorSpecification::SensorSpecification(
    const std::string &rType,
    const IndexType NewId)
    : IndexedObject(NewId),
      mType(rType),
      mSensorValue(0.0),
      mLocation{}
{
}

std::string SensorSpecification::GetType() const
{
    return mType;
}

void SensorSpecification::SetLocation(const array_1d<double, 3>& rLocation)
{
    mLocation = rLocation;
}

array_1d<double, 3> SensorSpecification::GetLocation() const
{
    return mLocation;
}

void SensorSpecification::SetSensorValue(const double SensorValue)
{
    mSensorValue = SensorValue;
}

double SensorSpecification::GetSensorValue() const
{
    return mSensorValue;
}

void SensorSpecification::AddNodalExpression(
    const std::string& rExpressionName,
    ContainerExpression<ModelPart::NodesContainerType>::Pointer pNodalExpression)
{
    KRATOS_TRY

    const auto p_itr = mNodalExpressions.find(rExpressionName);

    KRATOS_ERROR_IF_NOT(p_itr == mNodalExpressions.end())
        << "A nodal expression named \"" << rExpressionName << " already exists.";

    mNodalExpressions[rExpressionName] = pNodalExpression;

    KRATOS_CATCH("");
}

ContainerExpression<ModelPart::NodesContainerType>::Pointer SensorSpecification::GetNodalExpression(const std::string& rExpressionName) const
{
    KRATOS_TRY

    const auto p_itr = mNodalExpressions.find(rExpressionName);

    KRATOS_ERROR_IF(p_itr == mNodalExpressions.end())
        << "A nodal expression named \"" << rExpressionName << " not found.";

    return p_itr->second;

    KRATOS_CATCH("");
}

std::unordered_map<std::string, ContainerExpression<ModelPart::NodesContainerType>::Pointer> SensorSpecification::GetNodalExpressionsMap() const
{
    return mNodalExpressions;
}

void SensorSpecification::AddConditionExpression(
    const std::string& rExpressionName,
    ContainerExpression<ModelPart::ConditionsContainerType>::Pointer pConditionExpression)
{
    KRATOS_TRY

    const auto p_itr = mConditionExpressions.find(rExpressionName);

    KRATOS_ERROR_IF_NOT(p_itr == mConditionExpressions.end())
        << "A condition expression named \"" << rExpressionName << " already exists.";

    mConditionExpressions[rExpressionName] = pConditionExpression;

    KRATOS_CATCH("");
}

ContainerExpression<ModelPart::ConditionsContainerType>::Pointer SensorSpecification::GetConditionExpression(const std::string& rExpressionName) const
{
    KRATOS_TRY

    const auto p_itr = mConditionExpressions.find(rExpressionName);

    KRATOS_ERROR_IF(p_itr == mConditionExpressions.end())
        << "A condition expression named \"" << rExpressionName << " not found.";

    return p_itr->second;

    KRATOS_CATCH("");
}

std::unordered_map<std::string, ContainerExpression<ModelPart::ConditionsContainerType>::Pointer> SensorSpecification::GetConditionExpressionsMap() const
{
    return mConditionExpressions;
}

void SensorSpecification::AddElementExpression(
    const std::string& rExpressionName,
    ContainerExpression<ModelPart::ElementsContainerType>::Pointer pElementExpression)
{
    KRATOS_TRY

    const auto p_itr = mElementExpressions.find(rExpressionName);

    KRATOS_ERROR_IF_NOT(p_itr == mElementExpressions.end())
        << "A condition expression named \"" << rExpressionName << " already exists.";

    mElementExpressions[rExpressionName] = pElementExpression;

    KRATOS_CATCH("");
}

ContainerExpression<ModelPart::ElementsContainerType>::Pointer SensorSpecification::GetElementExpression(const std::string& rExpressionName) const
{
    KRATOS_TRY

    const auto p_itr = mElementExpressions.find(rExpressionName);

    KRATOS_ERROR_IF(p_itr == mElementExpressions.end())
        << "A condition expression named \"" << rExpressionName << " not found.";

    return p_itr->second;

    KRATOS_CATCH("");
}

std::unordered_map<std::string, ContainerExpression<ModelPart::ElementsContainerType>::Pointer> SensorSpecification::GetElementExpressionsMap() const
{
    return mElementExpressions;
}

std::vector<std::string> SensorSpecification::GetDataVariableNames() const
{
    std::vector<std::string> result;
    for (auto i_ptr = this->begin(); i_ptr != this->end(); ++i_ptr) {
        result.push_back(i_ptr->first->Name());
    }
    return result;
}

void SensorSpecification::ClearNodalExpressions()
{
    mNodalExpressions.clear();
}

void SensorSpecification::ClearConditionExpressions()
{
    mConditionExpressions.clear();
}

void SensorSpecification::ClearElementExpressions()
{
    mElementExpressions.clear();
}

std::string SensorSpecification::Info() const
{
    std::stringstream msg;
    msg << "SensorSpecification " << this->GetType() << " #" << this->GetId();
    return msg.str();
}

void SensorSpecification::PrintInfo(std::ostream& rOStream) const
{
    rOStream << Info();
}

void SensorSpecification::PrintData(std::ostream& rOStream) const
{
    PrintInfo(rOStream);
    rOStream << "    Value: " << this->GetSensorValue() << std::endl;
    rOStream << "    Location: " << this->GetLocation() << std::endl;
    DataValueContainer::PrintData(rOStream);
}

} /* namespace Kratos.*/