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
                *pContainer, *pVariable, [this, pVariable, StepIndex](auto& rValue, const Node& rNode) {
                    KRATOS_ERROR_IF_NOT(rNode.SolutionStepsDataHas(*pVariable))
                        << "The " << pVariable->Name() << " is not in the solution step variables list of "
                        << rNode << ".\n";

                    KRATOS_ERROR_IF_NOT(static_cast<int>(rNode.GetBufferSize()) >= this->mStepIndex)
                        << "The step index is larger than the nodal buffer size [ node buffer size = "
                        << rNode.GetBufferSize() << ", step index = " << this->mStepIndex
                        << ", node = " << rNode << " ].\n";

                    rValue = rNode.FastGetSolutionStepValue(*pVariable, StepIndex);
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

HistoricalVariableTensorAdaptor::BaseType::Pointer HistoricalVariableTensorAdaptor::Clone() const
{
    const auto& r_data_shape = this->DataShape();
    auto p_tensor_adaptor = Kratos::make_intrusive<HistoricalVariableTensorAdaptor>(mpContainer, mpVariable, std::vector<unsigned int>(r_data_shape.begin(), r_data_shape.end()));
    IndexPartition<IndexType>(p_tensor_adaptor->Size()).for_each([p_tensor_adaptor, this](const auto Index) {
        p_tensor_adaptor->ViewData()[Index] = this->ViewData()[Index];
    });
    return p_tensor_adaptor;
}

void HistoricalVariableTensorAdaptor::CollectData()
{
    std::visit(
        [this](auto pVariable) {
            TensorAdaptorUtils::CollectVariableData(
                this->Shape(), *this->mpContainer, *pVariable, this->ViewData(),
                [pVariable, this](auto& rValue, const Node& rNode) {
                    KRATOS_ERROR_IF_NOT(rNode.SolutionStepsDataHas(*pVariable))
                        << "The " << pVariable->Name() << " is not in the solution step variables list of "
                        << rNode << ".\n";

                    KRATOS_ERROR_IF_NOT(static_cast<int>(rNode.GetBufferSize()) >= this->mStepIndex)
                        << "The step index is larger than the nodal buffer size [ node buffer size = "
                        << rNode.GetBufferSize() << ", step index = " << this->mStepIndex
                        << ", node = " << rNode << " ].\n";

                    rValue = rNode.FastGetSolutionStepValue(*pVariable, this->mStepIndex);
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
                [pVariable, this](Node& rNode) -> auto& {
                    KRATOS_ERROR_IF_NOT(rNode.SolutionStepsDataHas(*pVariable))
                        << "The " << pVariable->Name() << " is not in the solution step variables list of "
                        << rNode << ".\n";

                    KRATOS_ERROR_IF_NOT(static_cast<int>(rNode.GetBufferSize()) >= this->mStepIndex)
                        << "The step index is larger than the nodal buffer size [ node buffer size = "
                        << rNode.GetBufferSize() << ", step index = " << this->mStepIndex
                        << ", node = " << rNode << " ].\n";

                    return rNode.FastGetSolutionStepValue(*pVariable, this->mStepIndex);
                });
        },
        mpVariable);
}

HistoricalVariableTensorAdaptor::ContainerPointerType HistoricalVariableTensorAdaptor::GetContainer() const
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