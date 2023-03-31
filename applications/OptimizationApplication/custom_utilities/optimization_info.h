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

    /**
     * @brief Advances the cyclic buffer by one index.
     *
     * The current instance and all the sub instances are advanced by one index
     * in the cyclic buffer.
     *
     */
    void AdvanceStep();

    /**
     * @brief Set the Buffer Size of the current Optimization info
     *
     * This sets the buffer size of the current optimization info instance. If @ref ResizeSubItems
     * is true, then all the sub items will also be resized to the new buffer size.
     *
     * @param BufferSize        New buffer size.
     * @param ResizeSubItems    If true, all sub items will also be resized.
     */
    void SetBufferSize(
        const IndexType BufferSize,
        const bool ResizeSubItems = false);

    /**
     * @brief Get the Buffer Size of the current instance
     *
     * @return IndexType        Buffer size
     */
    IndexType GetBufferSize() const;

    /**
     * @brief Checks whether given path is existing.
     *
     * This checks whether the given @ref rName exists in current instance or sub items.
     * The @ref rName should always starts with "/". Path is seperated by "/" character.
     *
     * When the lead name is found, then it checks whether the leaf name is present in the
     * respective container for the provided @ref StepIndex.
     *
     * An error is thrown if the leaf instance does not have the required buffer size.
     *
     * @param rName
     * @param StepIndex
     * @return true
     * @return false
     */
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

    IndexType GetBufferIndex(const IndexType StepIndex) const;

    ///@}
};

///@}
}