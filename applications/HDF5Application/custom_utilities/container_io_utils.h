//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  license: HDF5Application/license.txt
//
//  Main author:    Suneth Warnakulasuriya
//

// System includes

// Project includes
#include "containers/flags.h"
#include "includes/model_part.h"

// Application includes
#include "custom_utilities/data_type_utilities.h"

namespace Kratos
{
namespace HDF5
{
namespace Internals
{

class FlagIO
{
public:
    ///@name Public operations
    ///@{

    template<class TEntityType>
    inline char GetValue(
        const TEntityType& rEntity,
        const Flags& rFlag) const
    {
        return rEntity.Is(rFlag);
    }

    template<class TEntityType>
    inline void SetValue(
        TEntityType& rEntity,
        const Flags& rFlag,
        const char rValue) const
    {
        rEntity.Set(rFlag, static_cast<bool>(rValue));
    }

    ///@}
};

class HistoricalIO
{
public:
    ///@name Life cycle
    ///@{

    HistoricalIO(const unsigned int StepIndex) : mStepIndex(StepIndex) {}

    ///@}
    ///@name Public operations
    ///@{

    template<class TDataType>
    inline const TDataType& GetValue(
        const ModelPart::NodeType& rNode,
        const Variable<TDataType>& rVariable) const
    {
        return rNode.FastGetSolutionStepValue(rVariable, mStepIndex);
    }

    template<class TDataType>
    inline void SetValue(
        ModelPart::NodeType& rNode,
        const Variable<TDataType>& rVariable,
        const TDataType& rValue) const
    {
        rNode.FastGetSolutionStepValue(rVariable, mStepIndex) = rValue;
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    const unsigned int mStepIndex;

    ///@}
};

class NonHistoricalIO
{
public:
    ///@name Public operations
    ///@{

    template<class TEntityType, class TDataType>
    inline const TDataType& GetValue(
        const TEntityType& rEntity,
        const Variable<TDataType>& rVariable) const
    {
        return rEntity.GetValue(rVariable);
    }

    template<class TEntityType, class TDataType>
    inline void SetValue(
        TEntityType& rEntity,
        const Variable<TDataType>& rVariable,
        const TDataType& rValue) const
    {
        rEntity.SetValue(rVariable, rValue);
    }

    ///@}
};

template<class TContainerIOType>
std::string GetContainerIOType()
{
    if constexpr(std::is_same_v<TContainerIOType, HistoricalIO>) {
        return "HISTORICAL";
    } else if constexpr(std::is_same_v<TContainerIOType, NonHistoricalIO>) {
        return "NONHISTORICAL";
    } else if constexpr(std::is_same_v<TContainerIOType, FlagIO>) {
        return "FLAGS";
    } else {
        static_assert(!std::is_same_v<TContainerIOType, TContainerIOType>, "Unsupported container io type.");
    }
}

template<class TContainerType>
std::string GetContainerType()
{
    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        return "NODES";
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
        return "CONDITIONS";
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {
        return "ELEMENTS";
    } else {
        static_assert(!std::is_same_v<TContainerType, TContainerType>, "Unsupported container type.");
    }
}

template<class TContainerType, class TComponentType, class TContainerDataIO, class TDataType>
void CopyToContiguousDataArray(
    const TContainerType& rContainer,
    const TComponentType& rComponent,
    const TContainerDataIO& rContainerDataIO,
    TDataType* pBegin)
{
    KRATOS_TRY

    using value_type_traits = DataTypeTraits<typename ComponentTraits<TComponentType>::ValueType>;

    static_assert(value_type_traits::IsContiguous, "Only contiguous data types are supported.");

    if (rContainer.empty()) {
        // do nothing if the container is empty.
        return;
    }

    // get the first item for sizing.
    const auto& initial_value = rContainerDataIO.GetValue(rContainer.front(), rComponent);

    // get the stride from the first element to support dynamic types.
    const auto stride = value_type_traits::Size(initial_value);

    IndexPartition<unsigned int>(rContainer.size()).for_each([&rContainer, &rComponent, &rContainerDataIO, pBegin, stride](const auto Index) {
        const auto& value = rContainerDataIO.GetValue(*(rContainer.begin() + Index), rComponent);
        TDataType const* p_value_begin = value_type_traits::GetContiguousData(value);
        auto p_array_start = pBegin + Index * stride;
        std::copy(p_value_begin, p_value_begin + stride, p_array_start);
    });

    KRATOS_CATCH("");
}

template<class TContainerType, class TComponentType, class TContainerDataIO, class TDataType>
void CopyFromContiguousDataArray(
    TContainerType& rContainer,
    const TComponentType& rComponent,
    const TContainerDataIO& rContainerDataIO,
    TDataType const* pBegin,
    const std::vector<unsigned int>& rShape)
{
    KRATOS_TRY

    using value_type_traits = DataTypeTraits<typename ComponentTraits<TComponentType>::ValueType>;

    static_assert(value_type_traits::IsContiguous, "Only contiguous data types are supported.");

    using value_type = typename value_type_traits::ContainerType;

    if (rContainer.empty()) {
        // do nothing if the container is empty.
        return;
    }

    value_type tls_prototype;
    value_type_traits::Reshape(tls_prototype, rShape);

    const auto stride = value_type_traits::Size(tls_prototype);

    IndexPartition<unsigned int>(rContainer.size()).for_each(tls_prototype, [&rContainer, &rComponent, &rContainerDataIO, pBegin, stride](const auto Index, auto& rTLS) {
        TDataType * p_value_begin = value_type_traits::GetContiguousData(rTLS);
        TDataType const * p_array_start = pBegin + Index * stride;
        std::copy(p_array_start, p_array_start + stride, p_value_begin);
        rContainerDataIO.SetValue(*(rContainer.begin() + Index), rComponent, rTLS);
    });

    KRATOS_CATCH("");
}

} // namespace Internals
} // namespace HDF5
} // namespace Kratos