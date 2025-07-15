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

// External includes

// Project includes

// Include base h
#include "variable_tensor_adaptor.h"

namespace Kratos {

VariableTensorAdaptor::VariableTensorAdaptor(
    ContainerPointerType pContainer,
    VariablePointerType pVariable)
    : mpContainer(pContainer),
      mpVariable(pVariable)
{
    std::visit(
        [this](auto pContainer, auto pVariable) {
            this->SetShape(TensorAdaptorUtils::GetTensorShape(
                *pContainer, *pVariable, [pVariable](auto& rValue, const auto& rEntity) {
                    rValue = rEntity.GetValue(*pVariable);
                }));
        },
        mpContainer, mpVariable);
}

VariableTensorAdaptor::VariableTensorAdaptor(
    ContainerPointerType pContainer,
    VariablePointerType pVariable,
    const std::vector<unsigned int>& rDataShape)
    : mpContainer(pContainer),
      mpVariable(pVariable)
{
    std::visit(
        [&rDataShape, this](auto pContainer, auto pVariable) {
            this->SetShape(TensorAdaptorUtils::GetTensorShape(
                *pContainer, *pVariable, rDataShape.data(),
                rDataShape.data() + rDataShape.size()));
        },
        mpContainer, mpVariable);
}

VariableTensorAdaptor::BaseType::Pointer VariableTensorAdaptor::Clone() const
{
    return std::visit([this](auto pContainer) {
        const auto& r_data_shape = this->DataShape();
        auto p_tensor_adaptor = Kratos::make_intrusive<VariableTensorAdaptor>(pContainer, this->mpVariable, std::vector<unsigned int>(r_data_shape.begin(), r_data_shape.end()));
        IndexPartition<IndexType>(p_tensor_adaptor->Size()).for_each([p_tensor_adaptor, this](const auto Index) {
            p_tensor_adaptor->ViewData()[Index] = this->ViewData()[Index];
        });
        return p_tensor_adaptor;
    }, mpContainer);
}

void VariableTensorAdaptor::CollectData()
{
    std::visit(
        [this](auto pContainer, auto pVariable) {
            TensorAdaptorUtils::CheckAndSetValues(
                this->Shape(), *pContainer, *pVariable,
                [pVariable](const auto& rEntity) {
                    return rEntity.Has(*pVariable);
                },
                [pVariable](const auto& rValue, auto& rEntity) {
                    rEntity.SetValue(*pVariable, rValue);
                });

            TensorAdaptorUtils::CollectVariableData(
                this->Shape(), *pContainer, *pVariable, this->ViewData(),
                [pVariable](auto& rValue, const auto& rEntity) {
                    rValue = rEntity.GetValue(*pVariable);
                });
        },
        mpContainer, mpVariable);
}

void VariableTensorAdaptor::StoreData()
{
    std::visit(
        [this](auto pContainer, auto pVariable) {

            TensorAdaptorUtils::CheckAndSetValues(
                this->Shape(), *pContainer, *pVariable,
                [pVariable](const auto& rEntity) {
                    return rEntity.Has(*pVariable);
                },
                [pVariable](const auto& rValue, auto& rEntity) {
                    rEntity.SetValue(*pVariable, rValue);
                });

            TensorAdaptorUtils::StoreVariableData(
                this->Shape(), *pContainer, *pVariable, this->ViewData(),
                [pVariable](auto& rEntity) -> auto& {
                    return rEntity.GetValue(*pVariable);
                });
        },
        mpContainer, mpVariable);
}

VariableTensorAdaptor::ContainerPointerType VariableTensorAdaptor::GetContainer() const
{
    return mpContainer;
}

std::string VariableTensorAdaptor::Info() const
{
    return std::visit(
        [this](auto pContainer, auto pVariable) {
            return TensorAdaptorUtils::Info(
                "VariableTensorAdaptor: variable = " + pVariable->Name() + ", ",
                this->Shape(), *pContainer);
        },
        mpContainer, mpVariable);
}

} // namespace Kratos