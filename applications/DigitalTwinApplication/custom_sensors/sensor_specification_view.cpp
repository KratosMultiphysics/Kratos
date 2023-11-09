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
#include "sensor_specification_view.h"


namespace Kratos {

template<class TContainerType>
SensorSpecificationView<TContainerType>::SensorSpecificationView(
    SensorSpecification::Pointer pSensorSpecification,
    const std::string& rExpressionName)
    : mpSensorSpecification(pSensorSpecification)
{
    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        mpContainerExpression = mpSensorSpecification->GetNodalExpression(rExpressionName);
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
        mpContainerExpression = mpSensorSpecification->GetConditionExpression(rExpressionName);
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {
        mpContainerExpression = mpSensorSpecification->GetElementExpression(rExpressionName);
    } else {
        static_assert(std::is_same_v<TContainerType, TContainerType>, "Unsupported TContainerType.");
    }
}

template<class TContainerType>
SensorSpecification::Pointer SensorSpecificationView<TContainerType>::GetSensorSpecification() const
{
    return mpSensorSpecification;
}

template<class TContainerType>
typename ContainerExpression<TContainerType>::Pointer SensorSpecificationView<TContainerType>::GetContainerExpression() const
{
    return mpContainerExpression;
}

template<class TContainerType>
std::string SensorSpecificationView<TContainerType>::Info() const
{
    return mpSensorSpecification->Info();
}

template<class TContainerType>
void SensorSpecificationView<TContainerType>::PrintInfo(std::ostream& rOStream) const
{
    mpSensorSpecification->PrintInfo(rOStream);
}

template<class TContainerType>
void SensorSpecificationView<TContainerType>::PrintData(std::ostream& rOStream) const
{
    mpSensorSpecification->PrintData(rOStream);
}

// template insantiations
template class SensorSpecificationView<ModelPart::NodesContainerType>;
template class SensorSpecificationView<ModelPart::ConditionsContainerType>;
template class SensorSpecificationView<ModelPart::ElementsContainerType>;

} /* namespace Kratos.*/