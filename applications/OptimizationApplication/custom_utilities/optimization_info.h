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

    using IndexType = std::size_t;

    using OptimizationInfoType = OptimizationInfo<TArgs...>;

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
     *
     * When the leaf name is found, then it checks whether the leaf name is present in the
     * respective container for the provided @ref StepIndex.
     *
     * An error is thrown if the leaf instance does not have the required buffer size.
     *
     * @param rName                     Relative path to the value w.r.t. current instance.
     * @param StepIndex                 Step index of the cyclic buffer.
     * @return true                     if a value is existing.
     * @return false                    if a value is not existing.
     */
    bool HasValue(
        const std::string& rName,
        const IndexType StepIndex = 0) const;

    /**
     * @brief Checks the value is of the given type.
     *
     * This checks whether the value at @ref rName is of @ref TType.
     *
     * When the leaf name is found, then it checks whether the leaf name is present in the
     * respective container for the provided @ref StepIndex.
     *
     * An error is thrown if the leaf instance does not have the required buffer size.
     * An error is thrown if the path is not found.
     *
     * @tparam TType                    Type to be checked against in value.
     * @param rName                     Relative path to the value w.r.t. current instance.
     * @param StepIndex                 Step index of the cyclic buffer.
     * @return true                     If the value is of @ref TType.
     * @return false                    If the value is not of @ref TType.
     */
    template<class TType>
    bool IsValue(
        const std::string& rName,
        const IndexType StepIndex = 0) const;

    /**
     * @brief Get the value corresponding to rName.
     *
     * This returns the value corresponding to @ref rName path.
     *
     * An error is thrown if the leaf instance does not have the required buffer size.
     * An error is thrown if the path is not found.
     *
     * @param rName                     Relative path to the value w.r.t. current instance.
     * @param StepIndex                 Step index of the cyclic buffer.
     * @return ValueType                Return value.
     */
    ValueType GetValue(
        const std::string& rName,
        const IndexType StepIndex = 0) const;

    /**
     * @brief Get the reference value corresponding to rName.
     *
     * This returns the value reference corresponding to @ref rName path.
     *
     * An error is thrown if the leaf instance does not have the required buffer size.
     * An error is thrown if the path is not found.
     *
     * @param rName                     Relative path to the value w.r.t. current instance.
     * @param StepIndex                 Step index of the cyclic buffer.
     * @return ValueType                Return value.
     */
    ValueType& GetValue(
        const std::string& rName,
        const IndexType StepIndex = 0);

    /**
     * @brief Get the value corresponding to rName.
     *
     * This returns the value corresponding to @ref rName path.
     *
     * An error is thrown if the leaf instance does not have the required buffer size.
     * An error is thrown if the path is not found.
     * An error is thrown if the value found is not of the @ref TType.
     *
     * @tparam TType                    Type to be checked against in value.
     * @param rName                     Relative path to the value w.r.t. current instance.
     * @param StepIndex                 Step index of the cyclic buffer.
     * @return ValueType                Return value.
     */
    template<class TType>
    TType GetValue(
        const std::string& rName,
        const IndexType StepIndex = 0) const;

    /**
     * @brief Get the value reference corresponding to rName.
     *
     * This returns the value reference corresponding to @ref rName path.
     *
     * An error is thrown if the leaf instance does not have the required buffer size.
     * An error is thrown if the path is not found.
     * An error is thrown if the value found is not of the @ref TType.
     *
     * @tparam TType                    Type to be checked against in value.
     * @param rName                     Relative path to the value w.r.t. current instance.
     * @param StepIndex                 Step index of the cyclic buffer.
     * @return ValueType                Return value.
     */
    template<class TType>
    TType& GetValue(
        const std::string& rName,
        const IndexType StepIndex = 0);

    /**
     * @brief Set the Value at the given path
     *
     * This sets the given @ref rValue at the @ref rName.
     * Subitems will be created when reaching the leaf name from @ref rName.
     *
     * An error is thrown if a value already exists in the @ref rName path and Overwrite is false.
     * An error is thrown if the leaf instance does not have the required buffer size.
     *
     * @param rName                     Relative path to the value w.r.t. current instance.
     * @param rValue                    Value to be stored at rName path.
     * @param StepIndex                 Step index of the cyclic buffer.
     * @param Overwrite                 Whether to overwrite an existing value.
     */
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