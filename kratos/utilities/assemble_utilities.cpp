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
void AssembleUtilities::AssembleCurrentDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<int>& rVariable,
    const TMap<int>& rNodalValuesMap) const
{
    AssembleUtilities::CheckHistoricalVariable(rModelPart, rVariable);
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        rModelPart.Nodes(), rVariable, rNodalValuesMap,
        AssembleUtilities::UpdateHistoricalNodalValue<int>);
}

void AssembleUtilities::AssembleCurrentDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<double>& rVariable,
    const TMap<double>& rNodalValuesMap) const
{
    AssembleUtilities::CheckHistoricalVariable(rModelPart, rVariable);
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        rModelPart.Nodes(), rVariable, rNodalValuesMap,
        AssembleUtilities::UpdateHistoricalNodalValue<double>);
}

void AssembleUtilities::AssembleCurrentDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<array_1d<double, 3>>& rVariable,
    const TMap<array_1d<double, 3>>& rNodalValuesMap) const
{
    AssembleUtilities::CheckHistoricalVariable(rModelPart, rVariable);
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        rModelPart.Nodes(), rVariable, rNodalValuesMap,
        AssembleUtilities::UpdateHistoricalNodalValue<array_1d<double, 3>>);
}

void AssembleUtilities::AssembleNonHistoricalNodalDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<int>& rVariable,
    const TMap<int>& rNodalValuesMap) const
{
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        rModelPart.Nodes(), rVariable, rNodalValuesMap,
        AssembleUtilities::UpdateNonHistoricalValue<int, NodeType>);
}

void AssembleUtilities::AssembleNonHistoricalNodalDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<double>& rVariable,
    const TMap<double>& rNodalValuesMap) const
{
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        rModelPart.Nodes(), rVariable, rNodalValuesMap,
        AssembleUtilities::UpdateNonHistoricalValue<double, NodeType>);
}

void AssembleUtilities::AssembleNonHistoricalNodalDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<array_1d<double, 3>>& rVariable,
    const TMap<array_1d<double, 3>>& rNodalValuesMap) const
{
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        rModelPart.Nodes(), rVariable, rNodalValuesMap,
        AssembleUtilities::UpdateNonHistoricalValue<array_1d<double, 3>, NodeType>);
}

void AssembleUtilities::AssembleElementDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<int>& rVariable,
    const TMap<int>& rElementValuesMap) const
{
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        rModelPart.Elements(), rVariable, rElementValuesMap,
        AssembleUtilities::UpdateNonHistoricalValue<int, ElementType>);
}

void AssembleUtilities::AssembleElementDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<double>& rVariable,
    const TMap<double>& rElementValuesMap) const
{
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        rModelPart.Elements(), rVariable, rElementValuesMap,
        AssembleUtilities::UpdateNonHistoricalValue<double, ElementType>);
}

void AssembleUtilities::AssembleElementDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<array_1d<double, 3>>& rVariable,
    const TMap<array_1d<double, 3>>& rElementValuesMap) const
{
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        rModelPart.Elements(), rVariable, rElementValuesMap,
        AssembleUtilities::UpdateNonHistoricalValue<array_1d<double, 3>, ElementType>);
}

void AssembleUtilities::AssembleConditionDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<int>& rVariable,
    const TMap<int>& rConditionValuesMap) const
{
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        rModelPart.Conditions(), rVariable, rConditionValuesMap,
        AssembleUtilities::UpdateNonHistoricalValue<int, ConditionType>);
}

void AssembleUtilities::AssembleConditionDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<double>& rVariable,
    const TMap<double>& rConditionValuesMap) const
{
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        rModelPart.Conditions(), rVariable, rConditionValuesMap,
        AssembleUtilities::UpdateNonHistoricalValue<double, ConditionType>);
}

void AssembleUtilities::AssembleConditionDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<array_1d<double, 3>>& rVariable,
    const TMap<array_1d<double, 3>>& rConditionValuesMap) const
{
    AssembleUtilities::AssembleDataWithEntityValuesMap(
        rModelPart.Conditions(), rVariable, rConditionValuesMap,
        AssembleUtilities::UpdateNonHistoricalValue<array_1d<double, 3>, ConditionType>);
}

} // namespace Kratos