//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <vector>
#include <type_traits>

// External includes

// Project includes
#include "utilities/data_type_traits.h"
#include "utilities/container_io_utils.h"

// Include base h
#include "gauss_point_variable_tensor_adaptor.h"

namespace Kratos {

namespace GaussPointVariableTensorAdaptorHelperUtilities
{

template<class TContainerType, class TDataType>
DenseVector<unsigned int> GetTensorShape(
    const TContainerType& rContainer,
    const Variable<TDataType>& rVariable,
    const ProcessInfo& rProcessInfo)
{
    using entity_type = typename TContainerType::value_type;

    return TensorAdaptorUtils::GetTensorShape<std::vector<TDataType>>(
        rContainer, [&rProcessInfo, &rVariable](auto& rValues, const auto& rEntity) {
            const_cast<entity_type&>(rEntity).CalculateOnIntegrationPoints(
                rVariable, rValues, rProcessInfo);
        });
}

template<class TContainerType, class TDataType, class TSpanType, class TIntegerType>
void CollectData(
    const TContainerType& rContainer,
    const Variable<TDataType>& rVariable,
    const ProcessInfo& rProcessInfo,
    const TSpanType& rDataSpan,
    TIntegerType const * pShapeBegin,
    TIntegerType const * pShapeEnd)
{
    using entity_type = typename TContainerType::value_type;

    if constexpr(IsInList<TContainerType, ModelPart::ConditionsContainerType, ModelPart::ElementsContainerType>::value) {
        ContainerIOUtils::CopyToContiguousArray<std::vector<TDataType>>(
            rContainer, rDataSpan, pShapeBegin, pShapeEnd,
            [&rVariable, &rProcessInfo](auto& rValues, const auto& rEntity) {
                const_cast<entity_type&>(rEntity).CalculateOnIntegrationPoints(
                    rVariable, rValues, rProcessInfo);
            });
    }
}

template<class TContainerType, class TDataType, class TSpanType, class TIntegerType>
void StoreData(
    const TContainerType& rContainer,
    const Variable<TDataType>& rVariable,
    const ProcessInfo& rProcessInfo,
    const TSpanType& rDataSpan,
    TIntegerType const * pShapeBegin,
    TIntegerType const * pShapeEnd)
{
    if constexpr(IsInList<TContainerType, ModelPart::ConditionsContainerType, ModelPart::ElementsContainerType>::value) {
        ContainerIOUtils::CopyFromContiguousDataArray<std::vector<TDataType>>(
            rContainer, rDataSpan, pShapeBegin, pShapeEnd,
            [&rVariable, &rProcessInfo](const auto& rValues, auto& rEntity) {
                rEntity.SetValuesOnIntegrationPoints(rVariable, rValues, rProcessInfo);
            });
    }
}

} // namespace GaussPointVariableTensorAdaptorHelperUtilities

template<class TContainerPointerType>
GaussPointVariableTensorAdaptor::GaussPointVariableTensorAdaptor(
    TContainerPointerType pContainer,
    VariablePointerType pVariable,
    ProcessInfo::Pointer pProcessInfo)
    : mpContainer(pContainer),
      mpVariable(pVariable),
      mpProcessInfo(pProcessInfo)
{
    std::visit(
        [this, pContainer](auto pVariable) {
            this->SetShape(GaussPointVariableTensorAdaptorHelperUtilities::GetTensorShape(
                *pContainer, *pVariable, *this->mpProcessInfo));
        },
        mpVariable);
}

void GaussPointVariableTensorAdaptor::CollectData()
{
    std::visit(
        [this](auto pContainer, auto pVariable) {
            const auto& tensor_shape = this->Shape();
            GaussPointVariableTensorAdaptorHelperUtilities::CollectData(
                *pContainer, *pVariable, *this->mpProcessInfo, this->ViewData(),
                tensor_shape.data().begin(),
                tensor_shape.data().begin() + tensor_shape.size());
        },
        mpContainer, mpVariable);
}

void GaussPointVariableTensorAdaptor::StoreData()
{
    std::visit(
        [this](auto pContainer, auto pVariable) {
            const auto& tensor_shape = this->Shape();
            GaussPointVariableTensorAdaptorHelperUtilities::StoreData(
                *pContainer, *pVariable, *this->mpProcessInfo, this->ViewData(),
                tensor_shape.data().begin(),
                tensor_shape.data().begin() + tensor_shape.size());
        },
        mpContainer, mpVariable);
}

GaussPointVariableTensorAdaptor::ContainerPointerType GaussPointVariableTensorAdaptor::GetContainer() const
{
    return mpContainer;
}

std::string GaussPointVariableTensorAdaptor::Info() const
{
    return std::visit([this](auto pContainer, auto pVariable) {
        return TensorAdaptorUtils::Info("GaussPointVariableTensorAdaptor variable = " + pVariable->Name(), this->Shape(), *pContainer);
    }, mpContainer, mpVariable);
}

// template instantiations
template GaussPointVariableTensorAdaptor::GaussPointVariableTensorAdaptor(ModelPart::ConditionsContainerType::Pointer, VariablePointerType, ProcessInfo::Pointer);
template GaussPointVariableTensorAdaptor::GaussPointVariableTensorAdaptor(ModelPart::ElementsContainerType::Pointer, VariablePointerType, ProcessInfo::Pointer);

} // namespace Kratos