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

void VariableTensorAdaptor::CollectData()
{
    std::visit(
        [this](auto pContainer, auto pVariable) {
            TensorAdaptorUtils::CollectVariableData(
                this->Shape(), *pContainer, *pVariable, this->ViewData(),
                [pVariable](auto& rValue, const auto& rEntity) {
                    KRATOS_ERROR_IF_NOT(rEntity.Has(*pVariable)) << "The " << pVariable->Name() << " not found in the data value container of " << rEntity;
                    rValue = rEntity.GetValue(*pVariable);
                });
        },
        mpContainer, mpVariable);
}

void VariableTensorAdaptor::StoreData()
{
    std::visit(
        [this](auto pContainer, auto pVariable) {
            const auto& zero = TensorAdaptorUtils::GetZeroValue(*pVariable, this->DataShape());
            TensorAdaptorUtils::StoreVariableData(
                this->Shape(), *pContainer, *pVariable, this->ViewData(),
                [pVariable, &zero](auto& rEntity) -> auto& {
                    if (!rEntity.Has(*pVariable)) {
                        rEntity.SetValue(*pVariable, zero);
                    }
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