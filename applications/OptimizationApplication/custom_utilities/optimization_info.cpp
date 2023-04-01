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
        << GetBufferSize() << " [ StepIndex = " << StepIndex << " ]. OptimizationInfo:\n" << *this;

    return (mBufferIndex + StepIndex) % GetBufferSize();
}

template<class... TArgs>
bool OptimizationInfo<TArgs...>::HasValue(
    const std::string& rName,
    const IndexType StepIndex) const
{
    const auto delim_pos = rName.find('/');

    if (delim_pos == std::string::npos) {
        const auto& r_buffered_data = mBufferedData[GetBufferIndex(StepIndex)];
        return r_buffered_data.find(rName) != r_buffered_data.end() || mSubItems.find(rName) != mSubItems.end();
    } else {
        const auto sub_item_itr = mSubItems.find(rName.substr(0, delim_pos));
        if (sub_item_itr != mSubItems.end()) {
            return sub_item_itr->second->HasValue(rName.substr(delim_pos + 1), StepIndex);
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
    KRATOS_TRY

    const auto delim_pos = rName.find('/');

    if (delim_pos == std::string::npos) {
        if constexpr(std::is_same_v<TValueType, Pointer>) {
            const auto sub_item_itr = mSubItems.find(rName);
            return sub_item_itr != mSubItems.end();
        } else {
            const auto& r_buffered_data = mBufferedData[GetBufferIndex(StepIndex)];
            const auto data_itr = r_buffered_data.find(rName);
            if  (data_itr != r_buffered_data.end()) {
                return std::get_if<TValueType>(&data_itr->second) != nullptr;
            }
        }
    } else {
        const auto sub_item_itr = mSubItems.find(rName.substr(0, delim_pos));
        if (sub_item_itr != mSubItems.end()) {
            return sub_item_itr->second->template IsValue<TValueType>(rName.substr(delim_pos + 1), StepIndex);
        }
    }

    // throw an error if this block reaches this point, which means no value was returned.
    KRATOS_ERROR << "No value found for path \"" << rName << "\". OptimizationInfo:\n" << *this;

    return false;

    KRATOS_CATCH("");
}

template<class... TArgs>
template<class TValueType>
TValueType OptimizationInfo<TArgs...>::GetValue(
    const std::string& rName,
    const IndexType StepIndex) const
{
    KRATOS_TRY

    const auto delim_pos = rName.find('/');

    if (delim_pos == std::string::npos) {
        if constexpr(std::is_same_v<TValueType, Pointer>) {
            const auto sub_item_itr = mSubItems.find(rName);
            if (sub_item_itr != mSubItems.end()) {
                return sub_item_itr->second;
            }
        } else {
            const auto& r_buffered_data = mBufferedData[GetBufferIndex(StepIndex)];
            const auto data_itr = r_buffered_data.find(rName);
            if  (data_itr != r_buffered_data.end()) {
                const auto p_value = std::get_if<TValueType>(&data_itr->second);

                KRATOS_ERROR_IF(p_value == nullptr)
                    << "Found value at \"" << rName << "\" is not of the requested type. OptimizationInfo:\n" << *this;

                return *p_value;
            }
        }
    } else {
        const auto sub_item_itr = mSubItems.find(rName.substr(0, delim_pos));
        if (sub_item_itr != mSubItems.end()) {
            return sub_item_itr->second->template GetValue<TValueType>(rName.substr(delim_pos + 1), StepIndex);
        }
    }

    // throw an error if this block reaches this point, which means no value was returned.
    KRATOS_ERROR << "No value found for path \"" << rName << "\". OptimizationInfo:\n" << *this;

    return TValueType{};

    KRATOS_CATCH("");
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
        KRATOS_ERROR_IF_NOT(Overwrite || !HasValue(rName, StepIndex))
            << "A value already exists at \"" << rName << "\". OptimizationInfo:\n" << *this;

        auto& r_buffered_data = mBufferedData[GetBufferIndex(StepIndex)];
        auto& r_sub_items = this->mSubItems;

        std::visit([&r_buffered_data, &r_sub_items, &rName](const auto& rV) {
            OptimizationInfo<TArgs...>::AssignValue(r_buffered_data, r_sub_items, rName, rV);
        }, rValue);
    } else {
        const auto& r_name = rName.substr(0, delim_pos);
        auto sub_item_itr = mSubItems.find(r_name);
        if (sub_item_itr == mSubItems.end()) {
            KRATOS_ERROR_IF_NOT(Overwrite || !HasValue(r_name, StepIndex))
                << "A value already exists at \"" << r_name << "\". OptimizationInfo:\n" << *this;

            auto p_sub_item = std::make_shared<OptimizationInfoType>(this->GetBufferSize());
            mSubItems[r_name] = p_sub_item;
            p_sub_item->SetValue(rName.substr(delim_pos + 1), rValue, StepIndex, Overwrite);
        } else {
            sub_item_itr->second->SetValue(rName.substr(delim_pos + 1), rValue, StepIndex, Overwrite);
        }
    }

    KRATOS_CATCH("");
}

template<class... TArgs>
std::string OptimizationInfo<TArgs...>::Info(const std::string& rTab) const
{
    std::stringstream info;

    info << "{";

    for (IndexType i = 0; i < GetBufferSize(); ++i) {
        info << "\n" << rTab << "\t--- Step = " << i << " ---";
        for (const auto& r_buffer_item : mBufferedData[GetBufferIndex(i)]) {
            std::visit([&info, &rTab, &r_buffer_item](const auto& rValue) {
                info << "\n" << rTab << "\t\"" << r_buffer_item.first << "\": " << rValue << ",";
            }, r_buffer_item.second);
        }
    }

    std::stringstream tabbing;
    tabbing << rTab << "\t";
    for (const auto& r_sub_item : mSubItems) {
        info << "\n" << rTab << "\t\"" << r_sub_item.first << "\": " << r_sub_item.second->Info(tabbing.str()) << ",";
    }

    info << "\n" << rTab << "}";

    return info.str();
}

// template instantiation
template class OptimizationInfo<bool, int, double, std::string>;

template bool OptimizationInfo<bool, int, double, std::string>::IsValue<typename OptimizationInfo<bool, int, double, std::string>::Pointer>(const std::string&, const IndexType) const;
template typename OptimizationInfo<bool, int, double, std::string>::Pointer OptimizationInfo<bool, int, double, std::string>::GetValue<typename OptimizationInfo<bool, int, double, std::string>::Pointer>(const std::string&, const IndexType) const;

#define KRATOS_OPTIMIZATION_INFO_METHODS_FOR_TYPE(VALUE_TYPE)                                                                        \
    template bool OptimizationInfo<bool, int, double, std::string>::IsValue<VALUE_TYPE>(const std::string&, const IndexType) const;  \
    template VALUE_TYPE OptimizationInfo<bool, int, double, std::string>::GetValue<VALUE_TYPE>(const std::string&, const IndexType) const;

KRATOS_OPTIMIZATION_INFO_METHODS_FOR_TYPE(bool)
KRATOS_OPTIMIZATION_INFO_METHODS_FOR_TYPE(int)
KRATOS_OPTIMIZATION_INFO_METHODS_FOR_TYPE(double)
KRATOS_OPTIMIZATION_INFO_METHODS_FOR_TYPE(std::string)

#undef KRATOS_OPTIMIZATION_INFO_METHODS_FOR_TYPE

} // namespace Kratos