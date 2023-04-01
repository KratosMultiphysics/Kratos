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
#include "input_output/logger.h"
#include "utilities/string_utilities.h"

// Application includes

// Include base h
#include "optimization_info.h"

namespace Kratos {

template<class... TArgs>
OptimizationInfo<TArgs...>::OptimizationInfo(
    const IndexType BufferSize)
    : mBufferIndex(0),
      mBufferedData(),
      mSubItems()

{
    this->SetBufferSize(BufferSize);
}

template<class... TArgs>
void OptimizationInfo<TArgs...>::AdvanceStep()
{
    mBufferIndex = (mBufferIndex + 1) % GetBufferSize();

    for (auto& r_sub_item : mSubItems) {
        r_sub_item.second->AdvanceStep();
    }
}

template<class... TArgs>
void OptimizationInfo<TArgs...>::SetBufferSize(
    const IndexType BufferSize,
    const bool ResizeSubItems)
{
    if (mBufferedData.size() != BufferSize) {
        // first sets its own buffer
        mBufferedData.resize(BufferSize);

        // now sets the subitem buffers
        if (ResizeSubItems) {
            for (auto& r_sub_item : mSubItems) {
                r_sub_item.second->SetBufferSize(BufferSize);
            }
        }
    }
}

template<class... TArgs>
std::size_t OptimizationInfo<TArgs...>::GetBufferSize() const
{
    return mBufferedData.size();
}

template<class... TArgs>
std::size_t OptimizationInfo<TArgs...>::GetBufferIndex(const IndexType StepIndex) const
{
    KRATOS_ERROR_IF(StepIndex >= GetBufferSize())
        << "Invalid step index. Allowed step indices are < "
        << GetBufferSize() << " [ StepIndex = " << StepIndex << " ].\n";

    return (mBufferIndex - StepIndex) % GetBufferSize();
}

template<class... TArgs>
bool OptimizationInfo<TArgs...>::HasValue(
    const std::string& rName,
    const IndexType StepIndex) const
{
    const auto delim_pos = rName.find('/');

    if (delim_pos == std::string::npos) {
        const auto& r_buffered_data = mBufferedData[GetBufferIndex(StepIndex)];
        return r_buffered_data.find(rName) != r_buffered_data.end();
    } else {
        const auto sub_item_itr = mSubItems.find(rName.substr(0, delim_pos));
        if (sub_item_itr != mSubItems.end()) {
            return sub_item_itr->second->HasValue(rName.substr(delim_pos + 1));
        } else {
            return false;
        }
    }
}

template<class... TArgs>
template<class TValueType>
bool OptimizationInfo<TArgs...>::IsValue(
        const std::string& rName,
        const IndexType StepIndex) const
{
    // return GetValue(rName, StepIndex).has_value();
}

template<class... TArgs>
template<class TValueType>
std::optional<TValueType> OptimizationInfo<TArgs...>::GetValue(
    const std::string& rName,
    const IndexType StepIndex) const
{
    // KRATOS_TRY

    // const auto delim_pos = rName.find('/');

    // if (delim_pos == std::string::npos) {
    //     const auto& r_buffered_data = mBufferedData[GetBufferIndex(StepIndex)];
    //     const auto data_itr = r_buffered_data.find(rName);
    //     if  (data_itr != r_buffered_data.end()) {
    //         const auto p_value = std::get_if<TValueType>(data_itr->second);
    //         if (p_value != nullptr) {
    //             return *p_value;
    //         }
    //     } else {
    //         if constexpr(std::is_same_v<TValueType, OptimizationInfoPointer>) {
    //             const auto sub_item_itr = mSubItems.find(rName);
    //             if (sub_item_itr != mSubItems.end()) {
    //                 return sub_item_itr->second;
    //             }
    //         }
    //     }
    // } else {
    //     const auto sub_item_itr = mSubItems.find(rName.substr(0, delim_pos));
    //     if (sub_item_itr != mSubItems.end()) {
    //         return sub_item_itr->second->GetValue<TValueType>(rName.substr(delim_pos + 1));
    //     }
    // }

    // if (!is_found) {
    //     current_path << "\b";
    //     // put a nice error
    //     std::stringstream msg;
    //     msg << "The path \"" << current_path.str() << "\" not found. Parent path has following keys:";

    //     // first print the available buffered data
    //     if (StepIndex < p_optimization_info->mBufferedData.size()) {
    //         for (const auto& r_buffered_item : p_optimization_info->mBufferedData[StepIndex]) {
    //             msg <<"\n\t" << r_buffered_item.first;
    //         }
    //     }

    //     // now print the available sub_items
    //     for (const auto& r_sub_item : p_optimization_info->mSubItems) {
    //         msg << "\n\t" << r_sub_item.first;
    //     }

    //     KRATOS_ERROR << msg.str();
    // }

    // return std::nullopt;

    // KRATOS_CATCH("");
}

template<class... TArgs>
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
        auto& r_buffered_data = mBufferedData[GetBufferIndex(StepIndex)];
        std::visit([&r_buffered_data](const auto& rV) {
            // AssignValue(r_buffered_data, mSubItems, rName, rV);
        }, rValue);
    } else {
        const auto& r_name = rName.substr(0, delim_pos);
        auto sub_item_itr = mSubItems.find(r_name);
        if (sub_item_itr == mSubItems.end()) {
            auto p_sub_item = std::make_shared<OptimizationInfoType>(this->GetBufferSize());
            mSubItems[r_name] = p_sub_item;
            p_sub_item->SetValue(rName.substr(delim_pos + 1), rValue, StepIndex, Overwrite);
        } else {
            sub_item_itr->second->SetValue(rName.substr(delim_pos + 1), rValue, StepIndex, Overwrite);
        }
    }

    KRATOS_CATCH("");
}

// template instantiation
template class OptimizationInfo<bool, int, double, std::string>;

} // namespace Kratos