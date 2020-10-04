//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <unordered_map>
#include <vector>

// External includes

// Project includes
#include "containers/array_1d.h"
#include "containers/variable.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"

// Include base h
#include "utilities/assemble_utilities.h"

namespace Kratos
{
template <>
ModelPart::NodesContainerType& AssembleUtilities::GetContainer<ModelPart::NodesContainerType>(
    ModelPart& rModelPart)
{
    return rModelPart.Nodes();
}

template <>
ModelPart::ElementsContainerType& AssembleUtilities::GetContainer<ModelPart::ElementsContainerType>(
    ModelPart& rModelPart)
{
    return rModelPart.Elements();
}

template <>
ModelPart::ConditionsContainerType& AssembleUtilities::GetContainer<ModelPart::ConditionsContainerType>(
    ModelPart& rModelPart)
{
    return rModelPart.Conditions();
}

void AssembleUtilities::AssembleCurrentDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<int>& rVariable,
    const TMap<int>& rNodalValuesMap) const
{
    using entity_container_type = ModelPart::NodesContainerType;
    AssembleUtilities::CheckHistoricalVariable(rModelPart, rVariable);
    auto& container = AssembleUtilities::GetContainer<entity_container_type>(rModelPart);
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        container, rVariable, rNodalValuesMap,
        AssembleUtilities::UpdateHistoricalNodalValue<int>);
}

void AssembleUtilities::AssembleCurrentDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<double>& rVariable,
    const TMap<double>& rNodalValuesMap) const
{
    using entity_container_type = ModelPart::NodesContainerType;
    AssembleUtilities::CheckHistoricalVariable(rModelPart, rVariable);
    auto& container = AssembleUtilities::GetContainer<entity_container_type>(rModelPart);
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        container, rVariable, rNodalValuesMap,
        AssembleUtilities::UpdateHistoricalNodalValue<double>);
}

void AssembleUtilities::AssembleCurrentDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<array_1d<double, 3>>& rVariable,
    const TMap<array_1d<double, 3>>& rNodalValuesMap) const
{
    using entity_container_type = ModelPart::NodesContainerType;
    AssembleUtilities::CheckHistoricalVariable(rModelPart, rVariable);
    auto& container = AssembleUtilities::GetContainer<entity_container_type>(rModelPart);
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        container, rVariable, rNodalValuesMap,
        AssembleUtilities::UpdateHistoricalNodalValue<array_1d<double, 3>>);
}

void AssembleUtilities::AssembleNonHistoricalNodalDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<int>& rVariable,
    const TMap<int>& rNodalValuesMap) const
{
    using entity_container_type = ModelPart::NodesContainerType;
    auto& container = AssembleUtilities::GetContainer<entity_container_type>(rModelPart);
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        container, rVariable, rNodalValuesMap,
        AssembleUtilities::UpdateNonHistoricalValue<int, entity_container_type::value_type>);
}

void AssembleUtilities::AssembleNonHistoricalNodalDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<double>& rVariable,
    const TMap<double>& rNodalValuesMap) const
{
    using entity_container_type = ModelPart::NodesContainerType;
    auto& container = AssembleUtilities::GetContainer<entity_container_type>(rModelPart);
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        container, rVariable, rNodalValuesMap,
        AssembleUtilities::UpdateNonHistoricalValue<double, entity_container_type::value_type>);
}

void AssembleUtilities::AssembleNonHistoricalNodalDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<array_1d<double, 3>>& rVariable,
    const TMap<array_1d<double, 3>>& rNodalValuesMap) const
{
    using entity_container_type = ModelPart::NodesContainerType;
    auto& container = AssembleUtilities::GetContainer<entity_container_type>(rModelPart);
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        container, rVariable, rNodalValuesMap,
        AssembleUtilities::UpdateNonHistoricalValue<array_1d<double, 3>, entity_container_type::value_type>);
}

void AssembleUtilities::AssembleElementDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<int>& rVariable,
    const TMap<int>& rElementValuesMap) const
{
    using entity_container_type = ModelPart::ElementsContainerType;
    auto& container = AssembleUtilities::GetContainer<entity_container_type>(rModelPart);
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        container, rVariable, rElementValuesMap,
        AssembleUtilities::UpdateNonHistoricalValue<int, entity_container_type::value_type>);
}

void AssembleUtilities::AssembleElementDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<double>& rVariable,
    const TMap<double>& rElementValuesMap) const
{
    using entity_container_type = ModelPart::ElementsContainerType;
    auto& container = AssembleUtilities::GetContainer<entity_container_type>(rModelPart);
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        container, rVariable, rElementValuesMap,
        AssembleUtilities::UpdateNonHistoricalValue<double, entity_container_type::value_type>);
}

void AssembleUtilities::AssembleElementDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<array_1d<double, 3>>& rVariable,
    const TMap<array_1d<double, 3>>& rElementValuesMap) const
{
    using entity_container_type = ModelPart::ElementsContainerType;
    auto& container = AssembleUtilities::GetContainer<entity_container_type>(rModelPart);
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        container, rVariable, rElementValuesMap,
        AssembleUtilities::UpdateNonHistoricalValue<array_1d<double, 3>, entity_container_type::value_type>);
}

void AssembleUtilities::AssembleConditionDataWithValuesMap(ModelPart& rModelPart,
                                                           const Variable<int>& rVariable,
                                                           const TMap<int>& rConditionValuesMap) const
{
    using entity_container_type = ModelPart::ConditionsContainerType;
    auto& container = AssembleUtilities::GetContainer<entity_container_type>(rModelPart);
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        container, rVariable, rConditionValuesMap,
        AssembleUtilities::UpdateNonHistoricalValue<int, entity_container_type::value_type>);
}

void AssembleUtilities::AssembleConditionDataWithValuesMap(ModelPart& rModelPart,
                                                           const Variable<double>& rVariable,
                                                           const TMap<double>& rConditionValuesMap) const
{
    using entity_container_type = ModelPart::ConditionsContainerType;
    auto& container = AssembleUtilities::GetContainer<entity_container_type>(rModelPart);
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        container, rVariable, rConditionValuesMap,
        AssembleUtilities::UpdateNonHistoricalValue<double, entity_container_type::value_type>);
}

void AssembleUtilities::AssembleConditionDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<array_1d<double, 3>>& rVariable,
    const TMap<array_1d<double, 3>>& rConditionValuesMap) const
{
    using entity_container_type = ModelPart::ConditionsContainerType;
    auto& container = AssembleUtilities::GetContainer<entity_container_type>(rModelPart);
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        container, rVariable, rConditionValuesMap,
        AssembleUtilities::UpdateNonHistoricalValue<array_1d<double, 3>, entity_container_type::value_type>);
}

} // namespace Kratos