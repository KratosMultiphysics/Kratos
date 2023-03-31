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
        r_sub_item.second.AdvanceStep();
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
                r_sub_item.second.SetBufferSize(BufferSize);
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
    const auto& r_names = StringUtilities::SplitStringByDelimiter(rName, '/');

    OptimizationInfoType* p_optimization_info = this;
    for (IndexType i = 0; i < r_names.size(); ++i) {
        const auto& r_name = r_names[i];
        if (i == r_names.size() - 1) {
            // if the current index is the last, then it is the leaf
            auto& r_buffered_data = p_optimization_info->mBufferedData[GetBufferIndex(StepIndex)];
            return r_buffered_data.find(r_name) != r_buffered_data.end();
        } else {
            // it is not the last index, then this key should be present in the subitems
            auto sub_item_itr = p_optimization_info->mSubItems.find(r_name);
            if (sub_item_itr == p_optimization_info->mSubItems.end()) {
                return false;
            } else {
                p_optimization_info  = &sub_item_itr->second;
            }

        }
    }
}

template<class... TArgs>
template<class TType>
bool OptimizationInfo<TArgs...>::IsValue(
        const std::string& rName,
        const IndexType StepIndex) const
{
    return (std::get_if<TType>(this->GetValue(rName, StepIndex)) != nullptr);
}

template<class... TArgs>
typename OptimizationInfo<TArgs...>::ValueType OptimizationInfo<TArgs...>::GetValue(
    const std::string& rName,
    const IndexType StepIndex) const
{
    KRATOS_TRY

    const auto& r_names = StringUtilities::SplitStringByDelimiter(rName, '/');

    bool is_found = true;

    std::stringstream current_path;

    OptimizationInfoType* p_optimization_info = this;
    for (IndexType i = 0; i < r_names.size(); ++i) {
        const auto& r_name = r_names[i];
        current_path << r_name << "/";
        if (i == r_names.size() - 1) {
            // if the current index is the last, then it is the leaf
            auto& r_buffered_data = p_optimization_info->mBufferedData[GetBufferIndex(StepIndex)];
            auto sub_value = r_buffered_data.find(r_name);
            if (sub_value != r_buffered_data.end()) {
                return sub_value.second;
            } else {
                // now check whether this is a OptimizationInfo
                auto sub_item_itr = p_optimization_info->mSubItems.find(r_name);
                if (sub_item_itr != p_optimization_info->mSubItems.end()) {
                    return sub_item_itr.second;
                }
                is_found = false;
            }
        } else {
            // it is not the last index, then this key should be present in the subitems
            auto sub_item_itr = p_optimization_info->mSubItems.find(r_name);
            if (sub_item_itr == p_optimization_info->mSubItems.end()) {
                is_found = false;
                break;
            } else {
                p_optimization_info  = &sub_item_itr->second;
            }

        }
    }

    if (!is_found) {
        current_path << "\b";
        // put a nice error
        std::stringstream msg;
        msg << "The path \"" << current_path.str() << "\" not found. Parent path has following keys:";

        // first print the available buffered data
        if (StepIndex < p_optimization_info->mBufferedData.size()) {
            for (const auto& r_buffered_item : p_optimization_info->mBufferedData[StepIndex]) {
                msg <<"\n\t" r_buffered_item.first;
            }
        }

        // now print the available sub_items
        for (const auto& r_sub_item : p_optimization_info->mSubItems) {
            msg << "\n\t" << r_sub_item.first;
        }

        KRATOS_ERROR << msg.str();
    }

    KRATOS_CATCH("");
}

template<class... TArgs>
typename OptimizationInfo<TArgs...>::ValueType& OptimizationInfo<TArgs...>::GetValue(
    const std::string& rName,
    const IndexType StepIndex)
{
    KRATOS_TRY

    const auto& r_names = StringUtilities::SplitStringByDelimiter(rName, '/');

    bool is_found = true;

    std::stringstream current_path;

    OptimizationInfoType* p_optimization_info = this;
    for (IndexType i = 0; i < r_names.size(); ++i) {
        const auto& r_name = r_names[i];
        current_path << r_name << "/";
        if (i == r_names.size() - 1) {
            // if the current index is the last, then it is the leaf
            auto& r_buffered_data = p_optimization_info->mBufferedData[GetBufferIndex(StepIndex)];
            auto sub_value = r_buffered_data.find(r_name);
            if (sub_value != r_buffered_data.end()) {
                return sub_value.second;
            } else {
                // now check whether this is a OptimizationInfo
                auto sub_item_itr = p_optimization_info->mSubItems.find(r_name);
                if (sub_item_itr != p_optimization_info->mSubItems.end()) {
                    return sub_item_itr.second;
                }
                is_found = false;
            }
        } else {
            // it is not the last index, then this key should be present in the subitems
            auto sub_item_itr = p_optimization_info->mSubItems.find(r_name);
            if (sub_item_itr == p_optimization_info->mSubItems.end()) {
                is_found = false;
                break;
            } else {
                p_optimization_info  = &sub_item_itr->second;
            }

        }
    }

    if (!is_found) {
        current_path << "\b";
        // put a nice error
        std::stringstream msg;
        msg << "The path \"" << current_path.str() << "\" not found. Parent path has following keys:";

        // first print the available buffered data
        if (StepIndex < p_optimization_info->mBufferedData.size()) {
            for (const auto& r_buffered_item : p_optimization_info->mBufferedData[StepIndex]) {
                msg <<"\n\t" r_buffered_item.first;
            }
        }

        // now print the available sub_items
        for (const auto& r_sub_item : p_optimization_info->mSubItems) {
            msg << "\n\t" << r_sub_item.first;
        }

        KRATOS_ERROR << msg.str();
    }

    KRATOS_CATCH("");
}

template<class... TArgs>
template<class TType>
TType OptimizationInfo<TArgs...>::GetValue<TType>(
    const std::string& rName,
    const IndexType StepIndex ) const
{
    return std::get<TType>(GetValue(rName, StepIndex));
}

template<class... TArgs>
template<class TType>
TType& OptimizationInfo<TArgs...>::GetValue<TType>(
    const std::string& rName,
    const IndexType StepIndex)
{
    return std::get<TType>(GetValue(rName, StepIndex));
}

template<class... TArgs>
void OptimizationInfo<TArgs...>::SetValue(
    const std::string& rName,
    const ValueType& rValue,
    const IndexType StepIndex,
    const bool Overwrite)
{
    KRATOS_TRY

    const auto& r_names = StringUtilities::SplitStringByDelimiter(rName, '/');

    OptimizationInfoType* p_optimization_info = this;
    for (IndexType i = 0; i < r_names.size(); ++i) {
        const auto& r_name = r_names[i];
        if (i == r_names.size() - 1) {
            // if the current index is the last, then it is the leaf
            auto& r_buffered_data = p_optimization_info->mBufferedData[GetBufferIndex(StepIndex)];
            auto sub_value = r_buffered_data.find(r_name);
            KRATOS_ERROR_IF_NOT(Overwrite || sub_value != r_buffered_data.end()) << "A value at \"" << rName << "\" already exists.";
            r_buffered_data[r_name] = rValue;
        } else {
            // it is not the last index, then this key should be present in the subitems
            auto sub_item_itr = p_optimization_info->mSubItems.find(r_name);
            if (sub_item_itr == p_optimization_info->mSubItems.end()) {
                auto p_sub_item = std::make_shared<OptimizationInfoType>(this->GetBufferSize());
                p_optimization_info->mSubItems[r_name] = p_sub_item;
                p_optimization_info = p_sub_item.get();
            } else {
                p_optimization_info  = &sub_item_itr->second;
            }

        }
    }

    KRATOS_CATCH("");
}

// template instantiation
template class OptimizationInfo<bool, int, double, std::string>;

} // namespace Kratos