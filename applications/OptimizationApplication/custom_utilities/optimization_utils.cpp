//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

// System includes
#include <cmath>

// Project includes
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/variable_utils.h"

// Application includes

// Include base h
#include "optimization_utils.h"

namespace Kratos
{

template<class TContainerType>
GeometryData::KratosGeometryType OptimizationUtils::GetContainerEntityGeometryType(
    const TContainerType& rContainer,
    const DataCommunicator& rDataCommunicator)
{
    KRATOS_TRY

    if (rContainer.size() > 0) {
        const auto first_geometry_type = rContainer.begin()->GetGeometry().GetGeometryType();
        const bool local_value = block_for_each<MinReduction<bool>>(rContainer, [&](const auto& rEntity) -> IndexType {
            return first_geometry_type == rEntity.GetGeometry().GetGeometryType();
        });

        if (rDataCommunicator.AndReduceAll(local_value)) {
            return first_geometry_type;
        } else {
            return GeometryData::KratosGeometryType::Kratos_generic_type;
        }
    } else {
        return GeometryData::KratosGeometryType::Kratos_generic_type;
    }

    KRATOS_CATCH("")
}

template<class TContainerType, class TDataType>
bool OptimizationUtils::IsVariableExistsInAllContainerProperties(
    const TContainerType& rContainer,
    const Variable<TDataType>& rVariable,
    const DataCommunicator& rDataCommunicator)
{
    KRATOS_TRY

    const bool local_value = block_for_each<MinReduction<bool>>(rContainer, [&](const auto& rEntity){
        return rEntity.GetProperties().Has(rVariable);
    });

    return rDataCommunicator.AndReduceAll(local_value);

    KRATOS_CATCH("");
}

template<class TContainerType, class TDataType>
bool OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties(
    const TContainerType& rContainer,
    const Variable<TDataType>& rVariable,
    const DataCommunicator& rDataCommunicator)
{
    KRATOS_TRY

    const bool local_value = block_for_each<MaxReduction<bool>>(rContainer, [&](const auto& rEntity){
        return rEntity.GetProperties().Has(rVariable);
    });

    return rDataCommunicator.OrReduceAll(local_value);

    KRATOS_CATCH("");
}

template<class TContainerType>
void OptimizationUtils::CreateEntitySpecificPropertiesForContainer(
    ModelPart& rModelPart,
    TContainerType& rContainer)
{
    KRATOS_TRY

    auto element_properties_id = block_for_each<MaxReduction<IndexType>>(rContainer, [](const auto& rEntity) {
        return rEntity.GetProperties().Id();
    });

    auto properties_id = block_for_each<MaxReduction<IndexType>>(rModelPart.PropertiesArray(), [](const auto pProperties) {
        return pProperties->Id();
    });

    properties_id = std::max(element_properties_id, properties_id);

    // creation of properties is done in serial
    for (auto& r_entity : rContainer) {
        auto p_properties = rModelPart.CreateNewProperties(++properties_id);
        const auto& element_properties = r_entity.GetProperties();
        *p_properties = element_properties;
        p_properties->SetId(properties_id);
        r_entity.SetProperties(p_properties);
    }

    KRATOS_CATCH("");
}

template<>
KRATOS_API(OPTIMIZATION_APPLICATION) IndexType OptimizationUtils::GetVariableDimension(
    const Variable<double>& rVariable,
    const IndexType DomainSize)
{
    return 1;
}

template<>
KRATOS_API(OPTIMIZATION_APPLICATION) IndexType OptimizationUtils::GetVariableDimension(
    const Variable<array_1d<double, 3>>& rVariable,
    const IndexType DomainSize)
{
    return DomainSize;
}

void OptimizationUtils::CopySolutionStepVariablesList(
    ModelPart& rDestinationModelPart,
    const ModelPart& rOriginModelPart)
{
    rDestinationModelPart.GetNodalSolutionStepVariablesList() = rOriginModelPart.GetNodalSolutionStepVariablesList();
}

// template instantiations
template KRATOS_API(OPTIMIZATION_APPLICATION) GeometryData::KratosGeometryType OptimizationUtils::GetContainerEntityGeometryType(const ModelPart::ConditionsContainerType&, const DataCommunicator&);
template KRATOS_API(OPTIMIZATION_APPLICATION) GeometryData::KratosGeometryType OptimizationUtils::GetContainerEntityGeometryType(const ModelPart::ElementsContainerType&, const DataCommunicator&);

template KRATOS_API(OPTIMIZATION_APPLICATION) void OptimizationUtils::CreateEntitySpecificPropertiesForContainer(ModelPart&, ModelPart::ConditionsContainerType&);
template KRATOS_API(OPTIMIZATION_APPLICATION) void OptimizationUtils::CreateEntitySpecificPropertiesForContainer(ModelPart&, ModelPart::ElementsContainerType&);

template KRATOS_API(OPTIMIZATION_APPLICATION) bool OptimizationUtils::IsVariableExistsInAllContainerProperties(const ModelPart::ConditionsContainerType&, const Variable<double>&, const DataCommunicator& rDataCommunicator);
template KRATOS_API(OPTIMIZATION_APPLICATION) bool OptimizationUtils::IsVariableExistsInAllContainerProperties(const ModelPart::ElementsContainerType&, const Variable<double>&, const DataCommunicator& rDataCommunicator);
template KRATOS_API(OPTIMIZATION_APPLICATION) bool OptimizationUtils::IsVariableExistsInAllContainerProperties(const ModelPart::ConditionsContainerType&, const Variable<array_1d<double, 3>>&, const DataCommunicator& rDataCommunicator);
template KRATOS_API(OPTIMIZATION_APPLICATION) bool OptimizationUtils::IsVariableExistsInAllContainerProperties(const ModelPart::ElementsContainerType&, const Variable<array_1d<double, 3>>&, const DataCommunicator& rDataCommunicator);

template KRATOS_API(OPTIMIZATION_APPLICATION) bool OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties(const ModelPart::ConditionsContainerType&, const Variable<double>&, const DataCommunicator& rDataCommunicator);
template KRATOS_API(OPTIMIZATION_APPLICATION) bool OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties(const ModelPart::ElementsContainerType&, const Variable<double>&, const DataCommunicator& rDataCommunicator);
template KRATOS_API(OPTIMIZATION_APPLICATION) bool OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties(const ModelPart::ConditionsContainerType&, const Variable<array_1d<double, 3>>&, const DataCommunicator& rDataCommunicator);
template KRATOS_API(OPTIMIZATION_APPLICATION) bool OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties(const ModelPart::ElementsContainerType&, const Variable<array_1d<double, 3>>&, const DataCommunicator& rDataCommunicator);

}