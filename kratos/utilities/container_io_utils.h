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

#pragma once

// System includes
#include <vector>
#include <variant>
#include <sstream>
#include <type_traits>

// Project includes
#include "containers/flags.h"
#include "includes/ublas_interface.h"
#include "includes/model_part.h"
#include "utilities/data_type_traits.h"

// Application includes

namespace Kratos
{

class FlagIO
{
public:
    ///@name Type definitions
    ///@{

    using ReturnType = bool;

    ///@}
    ///@name Life cycle
    ///@{

    FlagIO(const Flags& rFlag) : mFlag(rFlag) {}

    ///@}
    ///@name Public operations
    ///@{

    template<class TEntityType>
    inline void GetValue(
        ReturnType& rOutput,
        const TEntityType& rEntity) const
    {
        rOutput = rEntity.Is(mFlag);
    }

    template<class TEntityType>
    inline void SetValue(
        const ReturnType& rInput,
        TEntityType& rEntity) const
    {
        rEntity.Set(mFlag, rInput);
    }

    std::string Info() const
    {
        return "Flag";
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    const Flags& mFlag;

    ///@}
};

template<class TDataType>
class HistoricalIO
{
public:
    ///@name Type definitions
    ///@{

    using DataType = TDataType;

    using ReturnType = TDataType;

    ///@}
    ///@name Life cycle
    ///@{

    HistoricalIO(
        const Variable<TDataType>& rVariable,
        const int StepIndex)
        : mpVariable(&rVariable),
          mStepIndex(StepIndex) {}

    ///@}
    ///@name Public operations
    ///@{

    inline void GetValue(
        ReturnType& rOutput,
        const Node& rNode) const
    {
        rOutput = rNode.FastGetSolutionStepValue(*mpVariable, mStepIndex);
    }

    inline void SetValue(
        const ReturnType& rInput,
        Node& rNode) const
    {
        rNode.FastGetSolutionStepValue(*mpVariable, mStepIndex) = rInput;
    }

    std::string Info() const
    {
        std::stringstream info;
        info << "Historical " << mpVariable->Name() << " at Step = " << mStepIndex;
        return info.str();
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    const Variable<TDataType>* mpVariable;

    const int mStepIndex;

    ///@}
};

template<class TDataType>
class NonHistoricalIO
{
public:
    ///@name Type definitions
    ///@{

    using DataType = TDataType;

    using ReturnType = TDataType;

    ///@}
    ///@name Life cycle
    ///@{

    NonHistoricalIO(const Variable<TDataType>& rVariable) : mpVariable(&rVariable) {}

    ///@}
    ///@name Public operations
    ///@{

    template<class TEntityType>
    inline void GetValue(
        ReturnType& rOutput,
        const TEntityType& rEntity) const
    {
        rOutput = rEntity.GetValue(*mpVariable);
    }

    template<class TEntityType>
    inline void SetValue(
        const ReturnType& rInput,
        TEntityType& rEntity) const
    {
        rEntity.SetValue(*mpVariable, rInput);
    }

    std::string Info() const
    {
        return mpVariable->Name();
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    const Variable<TDataType>* mpVariable;

    ///@}
};

template<class TDataType>
class GaussPointIO
{
public:
    ///@name Type definitions
    ///@{

    using DataType = TDataType;

    using ReturnType = std::vector<TDataType>;

    ///@}
    ///@name Life cycle
    ///@{

    GaussPointIO(
        const Variable<TDataType>& rVariable,
        const ProcessInfo& rProcessInfo)
        : mpVariable(&rVariable),
          mrProcessInfo(rProcessInfo) {}

    ///@}
    ///@name Public operations
    ///@{

    template<class TEntityType>
    inline void GetValue(
        ReturnType& rOutput,
        const TEntityType& rEntity) const
    {
        const_cast<TEntityType&>(rEntity).CalculateOnIntegrationPoints(*mpVariable, rOutput, mrProcessInfo);
    }

    template<class TEntityType>
    inline void SetValue(
        const ReturnType& rInput,
        TEntityType& rEntity) const
    {
        rEntity.SetValuesOnIntegrationPoints(*mpVariable, rInput, mrProcessInfo);
    }

    std::string Info() const
    {
        return "Gauss point " + mpVariable->Name();
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    const Variable<TDataType>* mpVariable;

    const ProcessInfo& mrProcessInfo;

    ///@}
};

template<class TContainerType, class TContainerDataIO>
void CopyToContiguousArray(
    const TContainerType& rContainer,
    const TContainerDataIO& rContainerDataIO,
    typename DataTypeTraits<typename TContainerDataIO::ReturnType>::PrimitiveType* pBegin,
    const IndexType Size)
{
    KRATOS_TRY

    using return_type = typename TContainerDataIO::ReturnType;

    using value_type_traits = DataTypeTraits<return_type>;

    if (rContainer.empty()) {
        // do nothing if the container is empty.
        return;
    }

    // get the first item for sizing.
    return_type initial_value;
    rContainerDataIO.GetValue(initial_value, rContainer.front());

    // get the stride from the first element to support dynamic types.
    const auto stride = value_type_traits::Size(initial_value);

    KRATOS_ERROR_IF_NOT(Size == rContainer.size() * stride)
        << "The contiguous array size mismatch with data in the container [ "
        << "Contiguous array size = " << Size << ", number of entities = "
        << rContainer.size() << ", data stride = " << stride << " ].";

    IndexPartition<unsigned int>(rContainer.size()).for_each(typename TContainerDataIO::ReturnType{}, [&rContainer,  &rContainerDataIO, stride, pBegin](const auto Index, auto& rTLS) {
        rContainerDataIO.GetValue(rTLS, *(rContainer.begin() + Index));
        auto p_subrange_begin = pBegin + Index * stride;
        value_type_traits::CopyToContiguousData(p_subrange_begin, rTLS);
    });

    KRATOS_CATCH("");
}

template<class TContainerType, class TContainerDataIO>
void CopyFromContiguousDataArray(
    TContainerType& rContainer,
    const TContainerDataIO& rContainerDataIO,
    typename DataTypeTraits<typename TContainerDataIO::ReturnType>::PrimitiveType const * pBegin,
    const std::vector<unsigned int>& rShape)
{
    KRATOS_TRY

    using return_type = typename TContainerDataIO::ReturnType;

    using value_type_traits = DataTypeTraits<return_type>;

    return_type dummy_value;
    value_type_traits::Reshape(dummy_value, rShape);

    const auto stride = value_type_traits::Size(dummy_value);

    IndexPartition<unsigned int>(rContainer.size()).for_each(dummy_value, [&rContainer, &rContainerDataIO, pBegin, stride](const auto Index, auto& rTLS) {
        auto p_subrange_begin = pBegin + Index * stride;
        value_type_traits::CopyFromContiguousData(rTLS, p_subrange_begin);
        rContainerDataIO.SetValue(rTLS, *(rContainer.begin() + Index));
    });

    KRATOS_CATCH("");
}

} // namespace Kratos