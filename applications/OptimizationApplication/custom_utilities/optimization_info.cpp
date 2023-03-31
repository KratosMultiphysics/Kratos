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
bool OptiimizationInfo<TArgs...>::DataItem::Has(const std::string& rName) const
{
    return mData.find(rName) != mData.end();
}

template<class... TArgs>
typename OptiimizationInfo<TArgs...>::ValueType OptiimizationInfo<TArgs...>::DataItem::GetValue(const std::string& rName) const
{
    return mData.find(rName)->second;
}

template<class... TArgs>
void OptiimizationInfo<TArgs...>::DataItem::SetValue(
    const std::string& rName,
    const ValueType& rValue)
{
    mData[rName] = rValue;
}

template<class... TArgs>
OptiimizationInfo<TArgs...>::OptiimizationInfo(const IndexType EchoLevel)
    : mEchoLevel(EchoLevel),
      mBufferSize(0),
      mBufferIndex(0),
      mRootDataItems()
{
}

template<class... TArgs>
void OptiimizationInfo<TArgs...>::SetBufferSize(const IndexType BufferSize)
{
    if (mRootDataItems.size() != BufferSize) {
        KRATOS_WARNING_IF("OptimizationInfo", mRootDataItems.size() != 0)
            << "Changing the buffer size from " << mRootDataItems.size() << " to "
            << BufferSize << " may lose the data in the current buffer.\n";

        mRootDataItems.resize(BufferSize);
        mBufferIndex = 0;

        KRATOS_INFO_IF("OptimizationInfo", mEchoLevel > 0)
            << "Resized optimization info buffer to " << BufferSize << ".\n";
    }
}

template<class... TArgs>
typename OptiimizationInfo<TArgs...>::ValueType OptiimizationInfo<TArgs...>::GetValue(
    const std::string& rName,
    const IndexType StepIndex) const
{
    KRATOS_TRY

    KRATOS_ERROR_IF(StepIndex >= mRootDataItems.size())
        << "Invalid step index. Allowed step indices are < "
        << mRootDataItems.size() << " [ StepIndex = " << StepIndex << " ].\n";

    const auto& r_names = StringUtilities::SplitStringByDelimiter(rName, '/');

    DataItem* p_data_item = &mRootDataItems[StepIndex];

    std::stringstream msg;
    for (const auto& r_name : r_names) {
        if (p_data_item->Has(r_name)) {



            std::visit([](auto Value) {
                if ()
            }, p_data_item->GetValue());
        }
    }


    KRATOS_CATCH("");
}

// template instantiation
template class OptiimizationInfo<bool, int, double, std::string>;

} // namespace Kratos