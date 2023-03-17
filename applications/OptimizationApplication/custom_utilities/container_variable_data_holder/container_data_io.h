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
#include <atomic>
#include <cmath>
#include <vector>

// Project includes
#include "containers/array_1d.h"
#include "includes/define.h"

// Application includes

namespace Kratos {

///@name Kratos Classes
///@{

class HistoricalContainerDataIO {
public:
    ///@name Public operations
    ///@{

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

    ///@}
};

class NonHistoricalContainerDataIO {
public:
    ///@name Public operations
    ///@{

    template <class TEntityType, class TDataType>
    static TDataType& GetValue(
        TEntityType& rEntity,
        const Variable<TDataType>& rVariable)
    {
        return rEntity.GetValue(rVariable);
    }

    template <class TEntityType, class TDataType>
    static void SetValue(
        TEntityType& rEntity,
        const Variable<TDataType>& rVariable,
        const TDataType& rValue)
    {
        rEntity.SetValue(rVariable, rValue);
    }

    ///@}
};

class PropertiesContainerDataIO {
public:
    ///@name Public operations
    ///@{

    template <class TEntityType, class TDataType>
    static TDataType& GetValue(
        TEntityType& rEntity,
        const Variable<TDataType>& rVariable)
    {
        return rEntity.GetProperties().GetValue(rVariable);
    }

    template <class TEntityType, class TDataType>
    static void SetValue(
        TEntityType& rEntity,
        const Variable<TDataType>& rVariable,
        const TDataType& rValue)
    {
        rEntity.GetProperties().SetValue(rVariable, rValue);
    }

    ///@}
};

} // namespace Kratos