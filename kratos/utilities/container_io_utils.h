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

// Project includes
#include "includes/define.h"
#include "containers/flags.h"
#include "includes/model_part.h"
#include "utilities/data_type_traits.h"

// Application includes

namespace Kratos
{

class FlagsIO
{
public:
    ///@name Type definitions
    ///@{

    using ReturnType = bool;

    template<class TContainerType>
    static constexpr bool IsAllowedContainer = IsInList<
                                                    TContainerType,
                                                    ModelPart::NodesContainerType,
                                                    ModelPart::ConditionsContainerType,
                                                    ModelPart::ElementsContainerType>::value;

    KRATOS_CLASS_POINTER_DEFINITION(FlagsIO);

    ///@}
    ///@name Life cycle
    ///@{

    FlagsIO(const Flags& rFlag) : mFlag(rFlag) {}

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

    using ReturnType = TDataType;

    template<class TContainerType>
    static constexpr bool IsAllowedContainer = IsInList<TContainerType, ModelPart::NodesContainerType>::value;

    KRATOS_CLASS_POINTER_DEFINITION(HistoricalIO);

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

    using ReturnType = TDataType;

    template<class TContainerType>
    static constexpr bool IsAllowedContainer = IsInList<TContainerType,
                                                        ModelPart::NodesContainerType,
                                                        ModelPart::ConditionsContainerType,
                                                        ModelPart::ElementsContainerType,
                                                        ModelPart::PropertiesContainerType,
                                                        ModelPart::GeometriesMapType,
                                                        ModelPart::MasterSlaveConstraintContainerType
                                                    >::value;

    KRATOS_CLASS_POINTER_DEFINITION(NonHistoricalIO);

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

    using ReturnType = std::vector<TDataType>;

    template<class TContainerType>
    static constexpr bool IsAllowedContainer = IsInList<TContainerType,
                                                        ModelPart::ConditionsContainerType,
                                                        ModelPart::ElementsContainerType
                                                        >::value;

    KRATOS_CLASS_POINTER_DEFINITION(GaussPointIO);

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

template<class TContainerType, class TContainerDataIO, class TIteratorType>
void CopyToContiguousArray(
    const TContainerType& rContainer,
    const TContainerDataIO& rContainerDataIO,
    typename DataTypeTraits<typename TContainerDataIO::ReturnType>::PrimitiveType* pBegin,
    TIteratorType pShapeBegin,
    TIteratorType pShapeEnd)
{
    KRATOS_TRY

    using return_type = typename TContainerDataIO::ReturnType;

    using value_type_traits = DataTypeTraits<return_type>;

    const auto stride = value_type_traits::Size(pShapeBegin + 1, pShapeEnd);

    IndexPartition<unsigned int>(rContainer.size()).for_each(typename TContainerDataIO::ReturnType{}, [&rContainer,  &rContainerDataIO, pShapeBegin, pShapeEnd, stride, pBegin](const auto Index, auto& rTLS) {
        rContainerDataIO.GetValue(rTLS, *(rContainer.begin() + Index));
        auto p_subrange_begin = pBegin + Index * stride;
        value_type_traits::CopyToContiguousData(p_subrange_begin, rTLS, pShapeBegin + 1, pShapeEnd);
    });

    KRATOS_CATCH("");
}

template<class TContainerType, class TContainerDataIO, class TIteratorType>
void CopyFromContiguousDataArray(
    TContainerType& rContainer,
    const TContainerDataIO& rContainerDataIO,
    typename DataTypeTraits<typename TContainerDataIO::ReturnType>::PrimitiveType const * pBegin,
    TIteratorType pShapeBegin,
    TIteratorType pShapeEnd)
{
    KRATOS_TRY

    using return_type = typename TContainerDataIO::ReturnType;

    using value_type_traits = DataTypeTraits<return_type>;

    return_type dummy_value{};
    if constexpr(DataTypeTraits<return_type>::Dimension > 0) {
        // skip for all the primitive types which does not need reshaping.
        value_type_traits::Reshape(dummy_value, &*(pShapeBegin + 1), &*(pShapeBegin) + std::distance(pShapeBegin, pShapeEnd));
    }

    const auto stride = value_type_traits::Size(pShapeBegin + 1, pShapeEnd);

    IndexPartition<unsigned int>(rContainer.size()).for_each(dummy_value, [&rContainer, &rContainerDataIO, pBegin, pShapeBegin, pShapeEnd, stride](const auto Index, auto& rTLS) {
        auto p_subrange_begin = pBegin + Index * stride;
        value_type_traits::CopyFromContiguousData(rTLS, p_subrange_begin, pShapeBegin + 1, pShapeEnd);
        rContainerDataIO.SetValue(rTLS, *(rContainer.begin() + Index));
    });

    KRATOS_CATCH("");
}

} // namespace Kratos