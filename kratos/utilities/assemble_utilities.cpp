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
template<>
ModelPart::NodesContainerType& AssembleUtilities::GetContainer<ModelPart::NodesContainerType>(
    ModelPart& rModelPart)
{
    return rModelPart.Nodes();
}

template<>
ModelPart::ElementsContainerType& AssembleUtilities::GetContainer<ModelPart::ElementsContainerType>(
    ModelPart& rModelPart)
{
    return rModelPart.Elements();
}

template<>
ModelPart::ConditionsContainerType& AssembleUtilities::GetContainer<ModelPart::ConditionsContainerType>(
    ModelPart& rModelPart)
{
    return rModelPart.Conditions();
}

void AssembleUtilities::AssembleCurrentDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<int>& rVariable,
    const TGPMap<NodeType, int>& rNodalValuesMap) const
{
    AssembleUtilities::CheckHistoricalVariable(rModelPart, rVariable);
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        rModelPart.Nodes(), rModelPart.GetCommunicator().GetDataCommunicator(), rVariable,
        rNodalValuesMap, AssembleUtilities::UpdateHistoricalNodalValue<int>);
    rModelPart.GetCommunicator().SynchronizeVariable(rVariable);
}

void AssembleUtilities::AssembleCurrentDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<double>& rVariable,
    const TGPMap<NodeType, double>& rNodalValuesMap) const
{
    AssembleUtilities::CheckHistoricalVariable(rModelPart, rVariable);
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        rModelPart.Nodes(), rModelPart.GetCommunicator().GetDataCommunicator(), rVariable,
        rNodalValuesMap, AssembleUtilities::UpdateHistoricalNodalValue<double>);
    rModelPart.GetCommunicator().SynchronizeVariable(rVariable);
}

void AssembleUtilities::AssembleCurrentDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<array_1d<double, 3>>& rVariable,
    const TGPMap<NodeType, array_1d<double, 3>>& rNodalValuesMap) const
{
    AssembleUtilities::CheckHistoricalVariable(rModelPart, rVariable);
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        rModelPart.Nodes(), rModelPart.GetCommunicator().GetDataCommunicator(),
        rVariable, rNodalValuesMap,
        AssembleUtilities::UpdateHistoricalNodalValue<array_1d<double, 3>>);
    rModelPart.GetCommunicator().SynchronizeVariable(rVariable);
}

void AssembleUtilities::AssembleNonHistoricalNodalDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<int>& rVariable,
    const TGPMap<NodeType, int>& rNodalValuesMap) const
{
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        rModelPart.Nodes(), rModelPart.GetCommunicator().GetDataCommunicator(), rVariable,
        rNodalValuesMap, AssembleUtilities::UpdateNonHistoricalValue<NodeType, int>);
    rModelPart.GetCommunicator().SynchronizeNonHistoricalVariable(rVariable);
}

void AssembleUtilities::AssembleNonHistoricalNodalDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<double>& rVariable,
    const TGPMap<NodeType, double>& rNodalValuesMap) const
{
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        rModelPart.Nodes(), rModelPart.GetCommunicator().GetDataCommunicator(),
        rVariable, rNodalValuesMap,
        AssembleUtilities::UpdateNonHistoricalValue<NodeType, double>);
    rModelPart.GetCommunicator().SynchronizeNonHistoricalVariable(rVariable);
}

void AssembleUtilities::AssembleNonHistoricalNodalDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<array_1d<double, 3>>& rVariable,
    const TGPMap<NodeType, array_1d<double, 3>>& rNodalValuesMap) const
{
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        rModelPart.Nodes(), rModelPart.GetCommunicator().GetDataCommunicator(),
        rVariable, rNodalValuesMap,
        AssembleUtilities::UpdateNonHistoricalValue<NodeType, array_1d<double, 3>>);
    rModelPart.GetCommunicator().SynchronizeNonHistoricalVariable(rVariable);
}

void AssembleUtilities::AssembleElementDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<int>& rVariable,
    const TGPMap<ElementType, int>& rElementValuesMap) const
{
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        rModelPart.Elements(), rModelPart.GetCommunicator().GetDataCommunicator(),
        rVariable, rElementValuesMap,
        AssembleUtilities::UpdateNonHistoricalValue<ElementType, int>);
}

void AssembleUtilities::AssembleElementDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<double>& rVariable,
    const TGPMap<ElementType, double>& rElementValuesMap) const
{
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        rModelPart.Elements(), rModelPart.GetCommunicator().GetDataCommunicator(),
        rVariable, rElementValuesMap,
        AssembleUtilities::UpdateNonHistoricalValue<ElementType, double>);
}

void AssembleUtilities::AssembleElementDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<array_1d<double, 3>>& rVariable,
    const TGPMap<ElementType, array_1d<double, 3>>& rElementValuesMap) const
{
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        rModelPart.Elements(), rModelPart.GetCommunicator().GetDataCommunicator(),
        rVariable, rElementValuesMap,
        AssembleUtilities::UpdateNonHistoricalValue<ElementType, array_1d<double, 3>>);
}

void AssembleUtilities::AssembleConditionDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<int>& rVariable,
    const TGPMap<ConditionType, int>& rConditionValuesMap) const
{
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        rModelPart.Conditions(), rModelPart.GetCommunicator().GetDataCommunicator(),
        rVariable, rConditionValuesMap,
        AssembleUtilities::UpdateNonHistoricalValue<ConditionType, int>);
}

void AssembleUtilities::AssembleConditionDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<double>& rVariable,
    const TGPMap<ConditionType, double>& rConditionValuesMap) const
{
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        rModelPart.Conditions(), rModelPart.GetCommunicator().GetDataCommunicator(),
        rVariable, rConditionValuesMap,
        AssembleUtilities::UpdateNonHistoricalValue<ConditionType, double>);
}

void AssembleUtilities::AssembleConditionDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<array_1d<double, 3>>& rVariable,
    const TGPMap<ConditionType, array_1d<double, 3>>& rConditionValuesMap) const
{
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        rModelPart.Conditions(), rModelPart.GetCommunicator().GetDataCommunicator(),
        rVariable, rConditionValuesMap,
        AssembleUtilities::UpdateNonHistoricalValue<ConditionType, array_1d<double, 3>>);
}

} // namespace Kratos