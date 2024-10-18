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
#include "expression/container_data_io.h"
#include "includes/model_part.h"

// Application includes

namespace Kratos {

///@name Kratos Classes
///@{

namespace ContainerDataIOTags {
    struct Properties    {};
} // namespace Tags

template <>
struct ContainerDataIO<ContainerDataIOTags::Properties>
{
    static constexpr std::string_view mInfo = "Properties";

    template<class TDataType, class TEntityType>
    static const TDataType& GetValue(
        const TEntityType& rEntity,
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