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
#include <string>
#include <variant>
#include <vector>
#include <unordered_map>

// Project includes
#include "includes/define.h"

// Application includes

namespace Kratos
{

///@name Kratos Classes
///@{

template<class... TArgs>
class KRATOS_API(OPTIMIZATION_APPLICATION) OptimizationInfo
{
public:
    ///@name Type definitions
    ///@{

    using OptimizationInfoType = OptimizationInfo<TArgs...>;

    using IndexType = std::size_t;

    using BufferedValueType = std::variant<TArgs...>;

    using BufferedMapType = std::unordered_map<std::string, BufferedValueType>;

    using BufferedDataType = std::vector<BufferedMapType>;

    using SubItemType = std::shared_ptr<OptimizationInfoType>;

    using SubItemsMap = std::unordered_map<std::string, SubItemType>;

    using ValueType = std::variant<SubItemType, TArgs...>;

    ///@}
    ///@name Life cycle
    ///@{

    OptimizationInfo(const IndexType BufferSize = 0);

    ///@}
    ///@name Public operations
    ///@{

    void SetBufferSize(
        const IndexType BufferSize,
        const bool ResizeSubItems = false);

    IndexType GetBufferSize() const;

    bool HasValue(
        const std::string& rName,
        const IndexType StepIndex = 0) const;

    template<class TType>
    bool IsValue(
        const std::string& rName,
        const IndexType StepIndex = 0) const;

    ValueType GetValue(
        const std::string& rName,
        const IndexType StepIndex = 0) const;

    void SetValue(
        const std::string& rName,
        const ValueType& rValue,
        const IndexType StepIndex = 0,
        const bool Overwrite = false);

    ///@}

private:
    ///@name Private member variables
    ///@{

    IndexType mBufferIndex;

    BufferedDataType mBufferedData;

    SubItemsMap mSubItems;

    ///@}
    ///@name Private operations
    ///@{

    void CheckStepIndex(const IndexType StepIndex) const;

    ///@}
};

///@}
}