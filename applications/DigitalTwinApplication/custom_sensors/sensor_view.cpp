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
#include <type_traits>

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
    : mpSensor(pSensor)
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
template class SensorView<ModelPart::NodesContainerType>;
template class SensorView<ModelPart::ConditionsContainerType>;
template class SensorView<ModelPart::ElementsContainerType>;

} /* namespace Kratos.*/