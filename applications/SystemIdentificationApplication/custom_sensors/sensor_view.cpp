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
    KRATOS_TRY

    auto p_expression = pSensor->GetContainerExpression(rExpressionName);
    this->mpContainerExpression = std::get<typename ContainerExpression<TContainerType>::Pointer>(p_expression);

    KRATOS_CATCH("");
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
    mpSensor->AddContainerExpression(name.str(), pContainerExpression);
}

template<class TContainerType>
typename ContainerExpression<TContainerType>::Pointer SensorView<TContainerType>::GetAuxiliaryExpression(const std::string& rSuffix) const
{
    std::stringstream name;
    name << this->mExpressionName << "_" << rSuffix;

    auto p_expression = mpSensor->GetContainerExpression(name.str());
    return std::get<typename ContainerExpression<TContainerType>::Pointer>(p_expression);
}

template<class TContainerType>
std::vector<std::string> SensorView<TContainerType>::GetAuxiliarySuffixes() const
{
    KRATOS_TRY

    std::vector<std::string> suffixes;
    for (const auto& r_pair : mpSensor->GetContainerExpressionsMap()) {
        const auto& r_name = r_pair.first;
        if (r_name.rfind(mExpressionName + "_", 0) == 0) {
            suffixes.push_back(r_name.substr(mExpressionName.size() + 1));
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

// template instantiations
template class SensorView<ModelPart::NodesContainerType>;
template class SensorView<ModelPart::ConditionsContainerType>;
template class SensorView<ModelPart::ElementsContainerType>;

} /* namespace Kratos.*/