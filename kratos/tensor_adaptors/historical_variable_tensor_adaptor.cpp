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
#include "historical_variable_tensor_adaptor.h"

namespace Kratos {

HistoricalVariableTensorAdaptor::HistoricalVariableTensorAdaptor(
    ModelPart::NodesContainerType::Pointer pContainer,
    VariablePointerType pVariable,
    const int StepIndex)
    : mpContainer(pContainer),
      mpVariable(pVariable),
      mStepIndex(StepIndex)
{
    std::visit(
        [this, pContainer, StepIndex](auto pVariable) {
            this->SetShape(TensorAdaptorUtils::GetTensorShape(
                *pContainer, *pVariable, [pVariable](auto& rValue, const auto& rEntity) {
                    rValue = rEntity.GetValue(*pVariable);
                }));
        },
        mpVariable);
}

HistoricalVariableTensorAdaptor::HistoricalVariableTensorAdaptor(
    ModelPart::NodesContainerType::Pointer pContainer,
    VariablePointerType pVariable,
    const std::vector<unsigned int>& rDataShape,
    const int StepIndex)
    : mpContainer(pContainer),
      mpVariable(pVariable),
      mStepIndex(StepIndex)
{
    std::visit([&rDataShape, this, pContainer](auto pVariable) {
        this->SetShape(TensorAdaptorUtils::GetTensorShape(
            *pContainer, *pVariable, rDataShape.data(),
            rDataShape.data() + rDataShape.size()));
    }, mpVariable);
}

void HistoricalVariableTensorAdaptor::CollectData()
{
    std::visit(
        [this](auto pVariable) {
            TensorAdaptorUtils::CollectVariableData(
                this->Shape(), *this->mpContainer, *pVariable, this->ViewData(),
                [pVariable, this](auto& rValue, const auto& rEntity) {
                    rValue = rEntity.FastGetSolutionStepValue(*pVariable, this->mStepIndex);
                });
        },
        mpVariable);
}

void HistoricalVariableTensorAdaptor::StoreData()
{
    std::visit(
        [this](auto pVariable) {
            TensorAdaptorUtils::StoreVariableData(
                this->Shape(), *this->mpContainer, *pVariable, this->ViewData(),
                [pVariable, this](const auto& rValue, auto& rEntity) {
                    rEntity.FastGetSolutionStepValue(*pVariable, this->mStepIndex) = rValue;
                });
        },
        mpVariable);
}

HistoricalVariableTensorAdaptor::ContainerType HistoricalVariableTensorAdaptor::GetContainer() const
{
    return mpContainer;
}

std::string HistoricalVariableTensorAdaptor::Info() const
{
    return std::visit(
        [this](auto pVariable) {
            return TensorAdaptorUtils::Info(
                "HistoricalVariableTensorAdaptor: variable = " + pVariable->Name() +
                    ", step index = " + std::to_string(this->mStepIndex) + ", ",
                this->Shape(), *this->mpContainer);
        },
        mpVariable);
}

} // namespace Kratos