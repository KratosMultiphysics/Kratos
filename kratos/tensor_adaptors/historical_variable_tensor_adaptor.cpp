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

namespace HistoricalVariableTensorAdaptorHelperUtils {

template<class TDataType>
void Check(
    const Node& rNode,
    const Variable<TDataType>& rVariable,
    const IndexType StepIndex)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(rNode.SolutionStepsDataHas(rVariable))
        << "The " << rVariable.Name() << " is not in the solution step variables list of "
        << rNode << ".\n";

    KRATOS_ERROR_IF_NOT(rNode.GetBufferSize() > StepIndex)
        << "The step index is larger than the nodal buffer size [ variable = "
        << rVariable.Name() << ", node buffer size = "
        << rNode.GetBufferSize() << ", step index = " << StepIndex
        << ", node = " << rNode << " ].\n";

    KRATOS_CATCH("");
}

} // namespace HistoricalVariableTensorAdaptorHelperUtils

HistoricalVariableTensorAdaptor::HistoricalVariableTensorAdaptor(
    ModelPart::NodesContainerType::Pointer pContainer,
    VariablePointerType pVariable,
    const int StepIndex)
    : mpVariable(pVariable),
      mStepIndex(StepIndex)
{
    std::visit([this, pContainer, StepIndex](auto pVariable) {
            this->mpStorage = Kratos::make_intrusive<TensorData<double>>(
                pContainer,
                TensorAdaptorUtils::GetTensorShape(
                    *pContainer, *pVariable, [pVariable, StepIndex](auto& rValue, const Node& rNode) {
                        HistoricalVariableTensorAdaptorHelperUtils::Check(rNode, *pVariable, StepIndex);
                        rValue = rNode.FastGetSolutionStepValue(*pVariable, StepIndex);
                    }));
        }, mpVariable);
}

HistoricalVariableTensorAdaptor::HistoricalVariableTensorAdaptor(
    ModelPart::NodesContainerType::Pointer pContainer,
    VariablePointerType pVariable,
    const std::vector<unsigned int>& rDataShape,
    const int StepIndex)
    : mpVariable(pVariable),
      mStepIndex(StepIndex)
{
    std::visit([&rDataShape, this, pContainer](auto pVariable) {
        this->mpStorage = Kratos::make_intrusive<TensorData<double>>(
                                pContainer, TensorAdaptorUtils::GetTensorShape(
                                *pContainer, *pVariable, rDataShape.data(),
                                rDataShape.data() + rDataShape.size()));
    }, mpVariable);
}

HistoricalVariableTensorAdaptor::HistoricalVariableTensorAdaptor(
    TensorData<double>::Pointer pTensorData,
    VariablePointerType pVariable,
    const int StepIndex)
    : mpVariable(pVariable),
      mStepIndex(StepIndex)
{
    this->mpStorage = pTensorData;

    // now check whether the given storage is compatible with the variable.
    std::visit([this](auto pVariable) {
        using data_type = BareType<decltype(*pVariable)>;
        const auto& r_data_shape = this->mpStorage->DataShape();
        KRATOS_ERROR_IF_NOT(DataTypeTraits<data_type>::IsValidShape(r_data_shape.data().begin(), r_data_shape.data().begin() + r_data_shape.size()))
            << "The data storage within the tensor data is not compatible with the " << pVariable->Name()
            << "[ tensor data = " << *(this->mpStorage) << " ].\n";

    }, mpVariable);

}

void HistoricalVariableTensorAdaptor::CollectData()
{
    std::visit([this](auto pContainer, auto pVariable) {
        using container_type = BareType<decltype(*pContainer)>;

        if constexpr(IsInList<container_type, ModelPart::NodesContainerType>) {
            using variable_type = BareType<decltype(*pVariable)>;
            using data_type = typename variable_type::Type;

            const auto& r_tensor_shape = this->Shape();

            ContainerIOUtils::CopyToContiguousArray<data_type>(
                *pContainer, this->ViewData(), r_tensor_shape.data().begin(),
                r_tensor_shape.data().begin() + r_tensor_shape.size(), [pVariable, this](auto& rValue, const Node& rNode) {
                        HistoricalVariableTensorAdaptorHelperUtils::Check(rNode, *pVariable, this->mStepIndex);
                        rValue = rNode.FastGetSolutionStepValue(*pVariable, this->mStepIndex);
                });
        }
    }, mpStorage->GetContainer(), mpVariable);
}

void HistoricalVariableTensorAdaptor::StoreData()
{
    std::visit([this](auto pContainer, auto pVariable){
        using container_type = BareType<decltype(*pContainer)>;

        if constexpr(IsInList<container_type, ModelPart::NodesContainerType>) {
            using variable_type = BareType<decltype(*pVariable)>;
            using data_type = typename variable_type::Type;

            const auto& r_tensor_shape = this->Shape();

            const auto& zero = TensorAdaptorUtils::GetZeroValue(*pVariable, this->DataShape());
            const auto& zero_shape = DataTypeTraits<data_type>::Shape(zero);

            ContainerIOUtils::CopyFromContiguousDataArray<data_type>(
                *pContainer, this->ViewData(), r_tensor_shape.data().begin(),
                r_tensor_shape.data().begin() + r_tensor_shape.size(), [&zero, &zero_shape, this, pVariable](Node& rNode) -> auto& {
                    HistoricalVariableTensorAdaptorHelperUtils::Check(rNode, *pVariable, this->mStepIndex);
                    auto& r_value = rNode.FastGetSolutionStepValue(*pVariable, this->mStepIndex);

                    // here we reshape the r_value to the given dimensions and sizes.
                    // The following method will not do anything if the type is static,
                    // but in the case where Variable<Vector> and Variable<Matrix>
                    // it will do the proper resizing.
                    // This adds no cost to the static data types.
                    if constexpr(DataTypeTraits<data_type>::IsDynamic) {
                        if (zero_shape != DataTypeTraits<data_type>::Shape(r_value)) {
                            // if the shape is not equal, then assign the correct shape.
                            r_value = zero;
                        }
                    }

                    return r_value;
                });
        }
    }, mpStorage->GetContainer(), mpVariable);
}

std::string HistoricalVariableTensorAdaptor::Info() const
{
    std::stringstream info;
    info << "HistoricalVariableTensorAdaptor:";
    std::visit([&info](auto pVariable) {
        info << " Variable = " << pVariable->Name();
    }, this->mpVariable);
    info << ", Step index = " << mStepIndex << ", " << *(this->mpStorage);
    return info.str();
}

} // namespace Kratos