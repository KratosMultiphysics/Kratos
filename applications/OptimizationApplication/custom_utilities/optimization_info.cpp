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

// System includes
#include <sstream>

// Project includes
#include "includes/model_part.h"
#include "containers/container_expression/container_data_io.h"
#include "containers/container_expression/specialized_container_expression.h"

// Application includes

// Include base h
#include "optimization_info.h"

namespace Kratos {

template <class... TArgs>
OptimizationInfo<TArgs...>::OptimizationInfo(const IndexType BufferSize)
    : mBufferIndex(0),
      mBufferedData(),
      mSubItems()

{
    this->SetBufferSize(BufferSize);
}

template <class... TArgs>
void OptimizationInfo<TArgs...>::AdvanceStep()
{
    // first advance the current instance
    mBufferIndex = (mBufferIndex + 1) % GetBufferSize();

    // now advance the sub instances
    for (auto& r_sub_item : mSubItems) {
        r_sub_item.second->AdvanceStep();
    }
}

template <class... TArgs>
void OptimizationInfo<TArgs...>::SetBufferSize(
    const IndexType BufferSize,
    const bool ResizeSubItems)
{
    // first sets its own buffer
    if (GetBufferSize() != BufferSize) {
        mBufferedData.resize(BufferSize);
    }

    // now sets the subitem buffers
    if (ResizeSubItems) {
        for (auto& r_sub_item : mSubItems) {
            r_sub_item.second->SetBufferSize(BufferSize, ResizeSubItems);
        }
    }
}

template <class... TArgs>
std::size_t OptimizationInfo<TArgs...>::GetBufferSize() const
{
    return mBufferedData.size();
}

template <class... TArgs>
std::size_t OptimizationInfo<TArgs...>::GetBufferIndex(const IndexType StepIndex) const
{
    KRATOS_ERROR_IF(StepIndex >= GetBufferSize())
        << "Invalid step index. Allowed step indices are < " << GetBufferSize()
        << " [ StepIndex = " << StepIndex << " ]. OptimizationInfo:\n"
        << *this;

    return (mBufferIndex + StepIndex) % GetBufferSize();
}

template <class... TArgs>
bool OptimizationInfo<TArgs...>::HasValue(
    const std::string& rName,
    const IndexType StepIndex) const
{
    const auto delim_pos = rName.find('/');

    if (delim_pos == std::string::npos) {
        // if the leaf name is found

        // get the buffer corresponding to StepIndex
        const auto& r_buffered_data = mBufferedData[GetBufferIndex(StepIndex)];

        // check the value in buffered data nd sub items
        return r_buffered_data.find(rName) != r_buffered_data.end() ||
               mSubItems.find(rName) != mSubItems.end();
    } else {
        // if not the leaf name, then check if it exists in sub items
        const auto sub_item_itr = mSubItems.find(rName.substr(0, delim_pos));
        if (sub_item_itr != mSubItems.end()) {
            // if sub item is found, call the same method recursively.
            return sub_item_itr->second->HasValue(rName.substr(delim_pos + 1), StepIndex);
        } else {
            // if the sub item is not found.
            return false;
        }
    }
}

template <class... TArgs>
template <class TValueType>
bool OptimizationInfo<TArgs...>::IsValue(
    const std::string& rName,
    const IndexType StepIndex) const
{
    KRATOS_TRY

    const auto delim_pos = rName.find('/');

    if (delim_pos == std::string::npos) {
        // if the leaf name is found
        if constexpr (std::is_same_v<TValueType, Pointer>) {
            // check if sub item is returned.
            const auto sub_item_itr = mSubItems.find(rName);
            // return only if the sub item is found. Otherwise, an error is thrown.
            if (sub_item_itr != mSubItems.end()) {
                return true;
            }
        } else {
            // if not sub item, then check the buffered data
            const auto& r_buffered_data = mBufferedData[GetBufferIndex(StepIndex)];
            const auto data_itr = r_buffered_data.find(rName);
            if (data_itr != r_buffered_data.end()) {
                // if buffered item found, check if it is the requested type.
                return std::get_if<TValueType>(&data_itr->second) != nullptr;
            }
        }
    } else {
        // if not the leaf name, then check if it exists in sub items
        const auto sub_item_itr = mSubItems.find(rName.substr(0, delim_pos));

        if (sub_item_itr != mSubItems.end()) {
            // if sub item exists, then call IsValue recursivly.
            return sub_item_itr->second->template IsValue<TValueType>(
                rName.substr(delim_pos + 1), StepIndex);
        }
    }

    if (HasValue(rName, StepIndex)) {
        // still no return is reached, then it is of a different type.
        return false;
    }

    // throw an error if this block reaches this point, which means no value was returned
    // and no value for given rName is found.
    KRATOS_ERROR << "No value found for path \"" << rName << "\". OptimizationInfo:\n"
                 << *this;

    return false;

    KRATOS_CATCH("");
}

template <class... TArgs>
template <class TValueType>
TValueType OptimizationInfo<TArgs...>::GetValue(
    const std::string& rName,
    const IndexType StepIndex) const
{
    KRATOS_TRY

    const auto delim_pos = rName.find('/');

    if (delim_pos == std::string::npos) {
        // if the leaf name is found
        if constexpr (std::is_same_v<TValueType, Pointer>) {
            // check if sub item is returned.
            const auto sub_item_itr = mSubItems.find(rName);
            // return only if the sub item is found. Otherwise, an error is thrown.
            if (sub_item_itr != mSubItems.end()) {
                return sub_item_itr->second;
            }
        } else {
            // if not sub item, then check the buffered data
            const auto& r_buffered_data = mBufferedData[GetBufferIndex(StepIndex)];
            const auto data_itr = r_buffered_data.find(rName);
            if (data_itr != r_buffered_data.end()) {
                // if buffered item found, check if it is the requested type.
                const auto p_value = std::get_if<TValueType>(&data_itr->second);

                KRATOS_ERROR_IF(p_value == nullptr)
                    << "Found value at \"" << rName
                    << "\" is not of the requested type. OptimizationInfo:\n"
                    << *this;

                return *p_value;
            }
        }
    } else {
        // if not the leaf name, then check if it exists in sub items
        const auto sub_item_itr = mSubItems.find(rName.substr(0, delim_pos));
        if (sub_item_itr != mSubItems.end()) {
            // if the sub item is found, then call GetValue iteratively.
            return sub_item_itr->second->template GetValue<TValueType>(
                rName.substr(delim_pos + 1), StepIndex);
        }
    }

    // if this is still reahced and value exists then.
    KRATOS_ERROR_IF(HasValue(rName, StepIndex))
        << "Found value at \"" << rName << "\" is not of the requested type. OptimizationInfo:\n"
        << *this;

    // throw an error if this block reaches this point, which means no value was returned
    // and no value for given rName is found.
    KRATOS_ERROR << "No value found for path \"" << rName << "\". OptimizationInfo:\n"
                 << *this;

    return TValueType{};

    KRATOS_CATCH("");
}

template <class... TArgs>
void OptimizationInfo<TArgs...>::SetValue(
    const std::string& rName,
    const ValueType& rValue,
    const IndexType StepIndex,
    const bool Overwrite)
{
    KRATOS_TRY

    const auto delim_pos = rName.find('/');

    // check if this is a leaf value
    if (delim_pos == std::string::npos) {
        // if the leaf name is found

        // first check whether Overwrite is specified and an existing value is found.
        // in both buffered data and sub items.
        KRATOS_ERROR_IF_NOT(Overwrite || !HasValue(rName, StepIndex))
            << "A value already exists at \"" << rName << "\". OptimizationInfo:\n"
            << *this;

        // get buffered data and sub item maps
        auto& r_buffered_data = mBufferedData[GetBufferIndex(StepIndex)];
        auto& r_sub_items = this->mSubItems;

        // assign to respective map the value
        std::visit(
            [&r_buffered_data, &r_sub_items, &rName, Overwrite](const auto& rV) {
                // first check whether alternate map has the same item,
                // and if so remove it from the alternate map.
                // It will be found in the Overwrite case only. Otherwise,
                // an error is thrown before.
                if (Overwrite) {
                    OptimizationInfo<TArgs...>::RemoveFromAlternateMap(
                        r_buffered_data, r_sub_items, rName, rV);
                }

                // now assign to respective map the given value.
                OptimizationInfo<TArgs...>::AssignValue(r_buffered_data,
                                                        r_sub_items, rName, rV);
            }, rValue);
    } else {
        // if rName is not a leaf
        const auto& r_name = rName.substr(0, delim_pos);
        auto sub_item_itr = mSubItems.find(r_name);
        if (sub_item_itr == mSubItems.end()) {
            // if a sub item is not found for the given name.
            // This check is required to see whether an item exists
            // in the buffered data. Otherwise, there can be
            // two items with same name in buffered data and
            // sub items.
            KRATOS_ERROR_IF_NOT(Overwrite || !HasValue(r_name, StepIndex))
                << "A value already exists at \"" << r_name << "\". OptimizationInfo:\n"
                << *this;

            // get buffered data and sub item maps
            auto& r_buffered_data = mBufferedData[GetBufferIndex(StepIndex)];
            auto& r_sub_items = this->mSubItems;

            if (Overwrite) {
                // first check whether alternate map has the same item,
                // and if so remove it from the alternate map.
                // It will be found in the Overwrite case only. Otherwise,
                // an error is thrown before.
                std::visit(
                    [&r_buffered_data, &r_sub_items, &r_name, Overwrite](const auto& rV) {
                        OptimizationInfo<TArgs...>::RemoveFromAlternateMap(
                            r_buffered_data, r_sub_items, r_name, rV);
                    }, rValue);
            }

            auto p_sub_item = std::make_shared<OptimizationInfoType>(this->GetBufferSize());
            mSubItems[r_name] = p_sub_item;
            p_sub_item->SetValue(rName.substr(delim_pos + 1), rValue, StepIndex, Overwrite);
        } else {
            sub_item_itr->second->SetValue(rName.substr(delim_pos + 1), rValue,
                                           StepIndex, Overwrite);
        }
    }

    KRATOS_CATCH("");
}

template <class... TArgs>
void OptimizationInfo<TArgs...>::GetKeys(
    std::vector<std::string>& rKeys,
    const IndexType StepIndex,
    const std::string& rPrefix) const
{
    // first put the buffered data
    if (StepIndex < GetBufferSize()) {
        const auto& r_buffered_data = mBufferedData[GetBufferIndex(StepIndex)];
        for (const auto& r_buffered_item_itr : r_buffered_data) {
            rKeys.push_back(std::string(rPrefix).append(r_buffered_item_itr.first));
        }
    }

    // now iterate through sub items and fill rKeys with leaves
    for (const auto& r_sub_item_itr : mSubItems) {
        r_sub_item_itr.second->GetKeys(
            rKeys, StepIndex,
            std::string(rPrefix).append(r_sub_item_itr.first).append("/"));
    }
}

template <class... TArgs>
std::string OptimizationInfo<TArgs...>::Info(const std::string& rTab) const
{
    std::stringstream info;

    info << "{";

    for (IndexType i = 0; i < GetBufferSize(); ++i) {
        info << "\n" << rTab << "\t--- Step = " << i << " ---";
        for (const auto& r_buffer_item : mBufferedData[GetBufferIndex(i)]) {
            std::visit(
                [&info, &rTab, &r_buffer_item](const auto& rValue) {
                    info << "\n"
                         << rTab << "\t\"" << r_buffer_item.first
                         << "\": " << rValue << ",";
                },
                r_buffer_item.second);
        }
    }

    std::stringstream tabbing;
    tabbing << rTab << "\t";
    for (const auto& r_sub_item : mSubItems) {
        info << "\n"
             << rTab << "\t\"" << r_sub_item.first
             << "\": " << r_sub_item.second->Info(tabbing.str()) << ",";
    }

    info << "\n" << rTab << "}";

    return info.str();
}

// template instantiation
#define KRATOS_OPTIMIZATION_INFO(...)                                                                                                                                                               \
    template class OptimizationInfo<__VA_ARGS__>;                                                                                                                                                   \
    template bool OptimizationInfo<__VA_ARGS__>::IsValue<typename OptimizationInfo<__VA_ARGS__>::Pointer>(const std::string&, const IndexType) const;                                               \
    template typename OptimizationInfo<__VA_ARGS__>::Pointer OptimizationInfo<__VA_ARGS__>::GetValue<typename OptimizationInfo<__VA_ARGS__>::Pointer>(const std::string&, const IndexType) const;

#define KRATOS_OPTIMIZATION_INFO_METHODS_FOR_BASIC_TYPE(VALUE_TYPE, ...)                                                                                                                            \
    template bool OptimizationInfo<__VA_ARGS__>::IsValue<VALUE_TYPE>( const std::string&, const IndexType) const;                                                                                   \
    template VALUE_TYPE OptimizationInfo<__VA_ARGS__>::GetValue<VALUE_TYPE>(const std::string&, const IndexType) const;

#define KRATOS_OPTIMIZATION_INFO_METHODS_FOR_CONTAINER_TYPE(CONTAINER_TYPE, CONTAINER_DATA_IO_TAG, ...)                                                                                                      \
    template bool OptimizationInfo<__VA_ARGS__>::IsValue<SpecializedContainerExpression<CONTAINER_TYPE, ContainerDataIO<CONTAINER_DATA_IO_TAG>>::Pointer>( const std::string&, const IndexType) const;       \
    template SpecializedContainerExpression<CONTAINER_TYPE, ContainerDataIO<CONTAINER_DATA_IO_TAG>>::Pointer OptimizationInfo<__VA_ARGS__>::GetValue<SpecializedContainerExpression<CONTAINER_TYPE, ContainerDataIO<CONTAINER_DATA_IO_TAG>>::Pointer>(const std::string&, const IndexType) const;

#define KRATOS_OPTIMIZATION_INFO_SUPPORTED_TYPES                                                                                         \
    ModelPart*,                                                                                                                          \
    bool,                                                                                                                                \
    int,                                                                                                                                 \
    double,                                                                                                                              \
    std::string,                                                                                                                         \
    SpecializedContainerExpression<ModelPart::NodesContainerType, ContainerDataIO<ContainerDataIOTags::Historical>>::Pointer,            \
    SpecializedContainerExpression<ModelPart::NodesContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>>::Pointer,         \
    SpecializedContainerExpression<ModelPart::ConditionsContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>>::Pointer,    \
    SpecializedContainerExpression<ModelPart::ElementsContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>>::Pointer

KRATOS_OPTIMIZATION_INFO(KRATOS_OPTIMIZATION_INFO_SUPPORTED_TYPES)

KRATOS_OPTIMIZATION_INFO_METHODS_FOR_BASIC_TYPE(bool, KRATOS_OPTIMIZATION_INFO_SUPPORTED_TYPES)
KRATOS_OPTIMIZATION_INFO_METHODS_FOR_BASIC_TYPE(int, KRATOS_OPTIMIZATION_INFO_SUPPORTED_TYPES)
KRATOS_OPTIMIZATION_INFO_METHODS_FOR_BASIC_TYPE(double, KRATOS_OPTIMIZATION_INFO_SUPPORTED_TYPES)
KRATOS_OPTIMIZATION_INFO_METHODS_FOR_BASIC_TYPE(std::string, KRATOS_OPTIMIZATION_INFO_SUPPORTED_TYPES)
KRATOS_OPTIMIZATION_INFO_METHODS_FOR_BASIC_TYPE(ModelPart*, KRATOS_OPTIMIZATION_INFO_SUPPORTED_TYPES)
KRATOS_OPTIMIZATION_INFO_METHODS_FOR_CONTAINER_TYPE(ModelPart::NodesContainerType, ContainerDataIOTags::Historical, KRATOS_OPTIMIZATION_INFO_SUPPORTED_TYPES)
KRATOS_OPTIMIZATION_INFO_METHODS_FOR_CONTAINER_TYPE(ModelPart::NodesContainerType, ContainerDataIOTags::NonHistorical, KRATOS_OPTIMIZATION_INFO_SUPPORTED_TYPES)
KRATOS_OPTIMIZATION_INFO_METHODS_FOR_CONTAINER_TYPE(ModelPart::ConditionsContainerType, ContainerDataIOTags::NonHistorical, KRATOS_OPTIMIZATION_INFO_SUPPORTED_TYPES)
KRATOS_OPTIMIZATION_INFO_METHODS_FOR_CONTAINER_TYPE(ModelPart::ElementsContainerType, ContainerDataIOTags::NonHistorical, KRATOS_OPTIMIZATION_INFO_SUPPORTED_TYPES)

#undef KRATOS_OPTIMIZATION_INFO_SUPPORTED_TYPES
#undef KRATOS_OPTIMIZATION_INFO_METHODS_FOR_BASIC_TYPE
#undef KRATOS_OPTIMIZATION_INFO_METHODS_FOR_CONTAINER_TYPE
#undef KRATOS_OPTIMIZATION_INFO

} // namespace Kratos