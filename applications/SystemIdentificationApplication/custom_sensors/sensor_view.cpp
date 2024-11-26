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
#include <sstream>

// External includes

// Project includes
#include "includes/model_part.h"

// Application includes

// Include base h
#include "sensor_view.h"


namespace Kratos {

template<class TContainerType>
SensorView<TContainerType>::SensorView(
    Sensor::Pointer pSensor,
    const std::string& rExpressionName)
    : mpSensor(pSensor),
      mExpressionName(rExpressionName)
{
    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        mpContainerExpression = mpSensor->GetNodalExpression(rExpressionName);
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
        mpContainerExpression = mpSensor->GetConditionExpression(rExpressionName);
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {
        mpContainerExpression = mpSensor->GetElementExpression(rExpressionName);
    } else {
        static_assert(std::is_same_v<TContainerType, TContainerType>, "Unsupported TContainerType.");
    }
}

template<class TContainerType>
Sensor::Pointer SensorView<TContainerType>::GetSensor() const
{
    return mpSensor;
}

template<class TContainerType>
typename ContainerExpression<TContainerType>::Pointer SensorView<TContainerType>::GetContainerExpression() const
{
    return mpContainerExpression;
}

template<class TContainerType>
std::string SensorView<TContainerType>::GetExpressionName() const
{
    return mExpressionName;
}

template<class TContainerType>
void  SensorView<TContainerType>::AddAuxiliaryExpression(
    const std::string& rSuffix,
    typename ContainerExpression<TContainerType>::Pointer pContainerExpression)
{
    std::stringstream name;
    name << this->mExpressionName << "_" << rSuffix;
    mpSensor->AddContainerExpression<TContainerType>(name.str(), pContainerExpression);
}

template<class TContainerType>
typename ContainerExpression<TContainerType>::Pointer SensorView<TContainerType>::GetAuxiliaryExpression(const std::string& rSuffix) const
{
    std::stringstream name;
    name << this->mExpressionName << "_" << rSuffix;

    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        return mpSensor->GetNodalExpression(name.str());
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
        return mpSensor->GetConditionExpression(name.str());
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {
        return mpSensor->GetElementExpression(name.str());
    } else {
        static_assert(std::is_same_v<TContainerType, TContainerType>, "Unsupported type.");
        return ContainerExpression<TContainerType>::Pointer;
    }
}

template<class TContainerType>
std::vector<std::string> SensorView<TContainerType>::GetAuxiliarySuffixes() const
{
    KRATOS_TRY

    std::vector<std::string> suffixes;
    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        for (const auto& r_pair : mpSensor->GetNodalExpressionsMap()) {
            const auto& r_name = r_pair.first;
            if (r_name.rfind(mExpressionName + "_", 0) == 0) {
                suffixes.push_back(r_name.substr(mExpressionName.size() + 1));
            }
        }
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
        for (const auto& r_pair : mpSensor->GetConditionExpressionsMap()) {
            const auto& r_name = r_pair.first;
            if (r_name.rfind(mExpressionName + "_", 0) == 0) {
                suffixes.push_back(r_name.substr(mExpressionName.size() + 1));
            }
        }
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {
        for (const auto& r_pair : mpSensor->GetElementExpressionsMap()) {
            const auto& r_name = r_pair.first;
            if (r_name.rfind(mExpressionName + "_", 0) == 0) {
                suffixes.push_back(r_name.substr(mExpressionName.size() + 1));
            }
        }
    }
    return suffixes;

    KRATOS_CATCH("");
}

template<class TContainerType>
std::string SensorView<TContainerType>::Info() const
{
    return mpSensor->Info();
}

template<class TContainerType>
void SensorView<TContainerType>::PrintInfo(std::ostream& rOStream) const
{
    mpSensor->PrintInfo(rOStream);
}

template<class TContainerType>
void SensorView<TContainerType>::PrintData(std::ostream& rOStream) const
{
    mpSensor->PrintData(rOStream);
}

// template insantiations
template KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) class SensorView<ModelPart::NodesContainerType>;
template KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) class SensorView<ModelPart::ConditionsContainerType>;
template KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) class SensorView<ModelPart::ElementsContainerType>;

} /* namespace Kratos.*/