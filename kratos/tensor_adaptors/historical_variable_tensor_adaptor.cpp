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
    const ModelPart::NodesContainerType& rContainer,
    const Variable<TDataType>& rVariable,
    const IndexType StepIndex)
{
    KRATOS_TRY

    // get the unique set of solution step variables lists
    std::set<VariablesList const *> solution_step_variables_lists;
    for (const auto& r_node : rContainer) {
        solution_step_variables_lists.insert(&*r_node.pGetVariablesList());

        KRATOS_ERROR_IF_NOT(r_node.GetBufferSize() > StepIndex)
            << "The step index is larger than the nodal buffer size [ variable = "
            << rVariable.Name() << ", node buffer size = "
            << r_node.GetBufferSize() << ", step index = " << StepIndex
            << ", node = " << r_node << " ].\n";
    }

    // now check whether the variable exists
    for (const auto& p_variables_list : solution_step_variables_lists) {
        KRATOS_ERROR_IF_NOT(p_variables_list->Has(rVariable))
            << "The " << rVariable.Name() << " is not in the solution step variables.\n";
    }

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
        HistoricalVariableTensorAdaptorHelperUtils::Check(*pContainer, *pVariable, StepIndex);
        this->mpStorage = Kratos::make_intrusive<Storage>(
            pContainer,
            TensorAdaptorUtils::GetTensorShape(
                *pContainer, *pVariable, [pVariable, StepIndex](auto& rValue, const Node& rNode) {
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
        this->mpStorage = Kratos::make_intrusive<Storage>(
                                pContainer, TensorAdaptorUtils::GetTensorShape(
                                *pContainer, *pVariable, rDataShape.data(),
                                rDataShape.data() + rDataShape.size()));
    }, mpVariable);
}

HistoricalVariableTensorAdaptor::HistoricalVariableTensorAdaptor(
    const TensorAdaptor& rOther,
    VariablePointerType pVariable,
    const int StepIndex,
    const bool Copy)
    : BaseType(rOther, Copy),
      mpVariable(pVariable),
      mStepIndex(StepIndex)
{
    // now check whether the given storage is compatible with the variable.
    std::visit([this](auto pVariable) {
        using data_type = typename BareType<decltype(*pVariable)>::Type;
        const auto& r_data_shape = this->mpStorage->DataShape();
        KRATOS_ERROR_IF_NOT(DataTypeTraits<data_type>::IsValidShape(r_data_shape.data().begin(), r_data_shape.data().begin() + r_data_shape.size()))
            << "The data storage within the tensor data is not compatible with the " << pVariable->Name()
            << "[ tensor data = " << this->mpStorage->Info() << " ].\n";

    }, mpVariable);

}

void HistoricalVariableTensorAdaptor::Check() const
{
    KRATOS_TRY

    std::visit([this](auto pContainer, auto pVariable) {
        using container_type = BareType<decltype(*pContainer)>;

        if constexpr(IsInList<container_type, ModelPart::NodesContainerType>) {
            const auto& r_tensor_shape = this->Shape();

            KRATOS_ERROR_IF_NOT(r_tensor_shape[0] == pContainer->size())
                << "Underlying container of the tensor data has changed size [ tensor data = "
                << this->mpStorage->Info() << ", container size = " << pContainer->size() << " ].\n";

            // first check if the variable is there, and step index is valid
            // This check is done every time CollectData or StoreData is called
            // because, the PointerVectorSet which the Storage holds
            // may have nodes from different model parts, or they may come from
            // a temporary PointerVectorSet which did not change in size, but
            // changed the underlying nodes or variables list.
            HistoricalVariableTensorAdaptorHelperUtils::Check(*pContainer, *pVariable, this->mStepIndex);
        }
    }, mpStorage->GetContainer(), mpVariable);

    KRATOS_CATCH("");
}

void HistoricalVariableTensorAdaptor::CollectData()
{

    std::visit([this](auto pContainer, auto pVariable) {
        using container_type = BareType<decltype(*pContainer)>;

        if constexpr(IsInList<container_type, ModelPart::NodesContainerType>) {
            using variable_type = BareType<decltype(*pVariable)>;
            using data_type = typename variable_type::Type;

            const auto& r_tensor_shape = this->Shape();

            KRATOS_ERROR_IF_NOT(r_tensor_shape[0] == pContainer->size())
                << "Underlying container of the tensor data has changed size [ tensor data = "
                << this->mpStorage->Info() << ", container size = " << pContainer->size() << " ].\n";

            ContainerIOUtils::CopyToContiguousArray<data_type>(
                *pContainer, this->ViewData(), r_tensor_shape.data().begin(),
                r_tensor_shape.data().begin() + r_tensor_shape.size(), [pVariable, this](auto& rValue, const Node& rNode) {
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

            KRATOS_ERROR_IF_NOT(r_tensor_shape[0] == pContainer->size())
                << "Underlying container of the tensor data has changed size [ tensor data = "
                << this->mpStorage->Info() << ", container size = " << pContainer->size() << " ].\n";

            if constexpr(DataTypeTraits<data_type>::IsDynamic) {
                const auto& zero = TensorAdaptorUtils::GetZeroValue(*pVariable, this->DataShape());
                const auto& zero_shape = DataTypeTraits<data_type>::Shape(zero);

                ContainerIOUtils::CopyFromContiguousDataArray<data_type>(
                    *pContainer, this->ViewData(), r_tensor_shape.data().begin(),
                    r_tensor_shape.data().begin() + r_tensor_shape.size(), [&zero, &zero_shape, this, pVariable](Node& rNode) -> auto& {
                        auto& r_value = rNode.FastGetSolutionStepValue(*pVariable, this->mStepIndex);

                        // here we reshape the r_value to the given dimensions and sizes.
                        // The following method will not do anything if the type is static,
                        // but in the case where Variable<Vector> and Variable<Matrix>
                        // it will do the proper resizing.
                        // This adds no cost to the static data types.
                        if (zero_shape != DataTypeTraits<data_type>::Shape(r_value)) {
                            // if the shape is not equal, then assign the correct shape.
                            r_value = zero;
                        }

                        return r_value;
                    });
            } else {
                ContainerIOUtils::CopyFromContiguousDataArray<data_type>(
                    *pContainer, this->ViewData(), r_tensor_shape.data().begin(),
                    r_tensor_shape.data().begin() + r_tensor_shape.size(), [this, pVariable](Node& rNode) -> auto& {
                        return rNode.FastGetSolutionStepValue(*pVariable, this->mStepIndex);
                    });
            }
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
    info << ", Step index = " << mStepIndex << ", " << this->mpStorage->Info();
    return info.str();
}

} // namespace Kratos