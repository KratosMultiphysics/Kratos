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

#pragma once

// System includes

// Project includes
#include "containers/variable.h"
#include "includes/model_part.h"

// Application includes

namespace Kratos {

///@name Kratos Classes
///@{

namespace ContainerDataIOTags {
    struct Historical    {};
    struct NonHistorical {};
    struct Properties    {};
} // namespace Tags

// Dummy generic class to be specialized later
template <class TContainerDataIOTags>
struct ContainerDataIO {};

template <>
struct ContainerDataIO<ContainerDataIOTags::Historical>
{
    template <class TDataType>
    static TDataType& GetValue(
        ModelPart::NodeType& rNode,
        const Variable<TDataType>& rVariable)
    {
        return rNode.FastGetSolutionStepValue(rVariable);
    }

    template <class TDataType>
    static void SetValue(
        ModelPart::NodeType& rNode,
        const Variable<TDataType>& rVariable,
        const TDataType& rValue)
    {
        rNode.FastGetSolutionStepValue(rVariable) = rValue;
    }
};

template <>
struct ContainerDataIO<ContainerDataIOTags::NonHistorical>
{
    template<class TDataType, class TEntityType>
    static TDataType& GetValue(
        TEntityType& rEntity,
        const Variable<TDataType>& rVariable)
    {
        return rEntity.GetValue(rVariable);
    }

    template<class TDataType, class TEntityType>
    static void SetValue(
        TEntityType& rEntity,
        const Variable<TDataType>& rVariable,
        const TDataType& rValue)
    {
        rEntity.SetValue(rVariable, rValue);
    }
};

template <>
struct ContainerDataIO<ContainerDataIOTags::Properties>
{
    template<class TDataType, class TEntityType>
    static TDataType& GetValue(
        TEntityType& rEntity,
        const Variable<TDataType>& rVariable)
    {
        static_assert(!(std::is_same_v<TEntityType, ModelPart::NodeType>), "Properties retrieval is only supported for element and conditions.");
        return rEntity.GetProperties().GetValue(rVariable);
    }

    template<class TDataType, class TEntityType>
    static void SetValue(
        TEntityType& rEntity,
        const Variable<TDataType>& rVariable,
        const TDataType& rValue)
    {
        static_assert(!(std::is_same_v<TEntityType, ModelPart::NodeType>), "Properties setter is only supported for element and conditions.");
        rEntity.GetProperties().SetValue(rVariable, rValue);
    }
};

} // namespace Kratos