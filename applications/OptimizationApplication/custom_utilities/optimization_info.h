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
// #include "inclu"

// Application includes

namespace Kratos
{

///@name Kratos Classes
///@{

template<class... TArgs>
class KRATOS_API(OPTIMIZATION_APPLICATION) OptiimizationInfo
{
private:
    ///@name Private class forward declarations
    ///@{

    class DataItem;

    ///@}
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using ValueType = std::variant<std::shared_ptr<DataItem>, TArgs...>;

    using MapType = std::unordered_map<std::string, ValueType>;

    ///@}
    ///@name Life cycle
    ///@{

    OptiimizationInfo(const IndexType EchoLevel = 0);

    ///@}
    ///@name Public operations
    ///@{

    void SetBufferSize(const IndexType BufferSize);

    ValueType GetValue(
        const std::string& rName,
        const IndexType StepIndex = 0) const;

    ///@}

private:
    ///@name Private class declaration
    ///@{

    class DataItem
    {
    public:
        bool Has(const std::string& rName) const;

        ValueType GetValue(
            const std::string& rName) const;

        void SetValue(
            const std::string& rName,
            const ValueType& rValue);

    private:
        MapType mData;
    };

    ///@}
    ///@name Private member variables
    ///@{

    const IndexType mEchoLevel;

    IndexType mBufferSize;

    IndexType mBufferIndex;

    std::vector<DataItem> mRootDataItems;

    ///@}
};

///@}
}