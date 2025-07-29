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
    : mpVariable(pVariable)
{
    std::visit([this](auto pContainer, auto pVariable) {
        this->mpStorage = Kratos::make_intrusive<TensorData<double>>(
                                pContainer, TensorAdaptorUtils::GetTensorShape(
                                                *pContainer, *pVariable,
                                                [pVariable](auto& rValue, const auto& rEntity) {
                                                    rValue = rEntity.GetValue(*pVariable);
                                                }));
    }, pContainer, mpVariable);
}

VariableTensorAdaptor::VariableTensorAdaptor(
    ContainerPointerType pContainer,
    VariablePointerType pVariable,
    const std::vector<unsigned int>& rDataShape)
    : mpVariable(pVariable)
{
    std::visit([&rDataShape, this](auto pContainer, auto pVariable) {
        this->mpStorage = Kratos::make_intrusive<TensorData<double>>(
                                pContainer, TensorAdaptorUtils::GetTensorShape(
                                *pContainer, *pVariable, rDataShape.data(),
                                rDataShape.data() + rDataShape.size()));
    }, pContainer, mpVariable);
}

VariableTensorAdaptor::VariableTensorAdaptor(
    TensorData<double>::Pointer pTensorData,
    VariablePointerType pVariable)
    : mpVariable(pVariable)
{
    KRATOS_TRY

    this->mpStorage = pTensorData;

    // now check whether the given storage is compatible with the variable.
    std::visit([this](auto pVariable) {
        using data_type = BareType<decltype(*pVariable)>;
        const auto& r_data_shape = this->mpStorage->DataShape();
        KRATOS_ERROR_IF_NOT(DataTypeTraits<data_type>::IsValidShape(r_data_shape.data().begin(), r_data_shape.data().begin() + r_data_shape.size()))
            << "The data storage within the tensor data is not compatible with the " << pVariable->Name()
            << "[ tensor data = " << *(this->mpStorage) << " ].\n";

    }, mpVariable);

    KRATOS_CATCH("");
}

void VariableTensorAdaptor::CollectData()
{
    std::visit([this](auto pContainer, auto pVariable) {
        using variable_type = BareType<decltype(*pVariable)>;
        using data_type = typename variable_type::Type;

        const auto& r_tensor_shape = this->Shape();

        KRATOS_ERROR_IF_NOT(r_tensor_shape[0] == pContainer->size())
            << "Underlying container of the tensor data has changed size [ tensor data = "
            << *this->GetTensorData() << ", container size = " << pContainer->size() << " ].\n";

        ContainerIOUtils::CopyToContiguousArray<data_type>(
            *pContainer, this->ViewData(), r_tensor_shape.data().begin(),
            r_tensor_shape.data().begin() + r_tensor_shape.size(),
            [pVariable](auto& rValue, const auto& rEntity) {
                auto p_value = rEntity.pGetValue(*pVariable);
                KRATOS_ERROR_IF_NOT(p_value)
                    << "The " << pVariable->Name()
                    << " not found in the data value container of " << rEntity;
                rValue = *p_value;
            });
    }, this->mpStorage->GetContainer(), mpVariable);
}

void VariableTensorAdaptor::StoreData()
{
    std::visit([this](auto pContainer, auto pVariable) {
        using variable_type = BareType<decltype(*pVariable)>;
        using data_type = typename variable_type::Type;

        const auto& r_tensor_shape = this->Shape();

        KRATOS_ERROR_IF_NOT(r_tensor_shape[0] == pContainer->size())
            << "Underlying container of the tensor data has changed size [ tensor data = "
            << *this->GetTensorData() << ", container size = " << pContainer->size() << " ].\n";

        // this zero value may be different from the Variable::Zero()
        // in the case where the variable type is Variable<Vector> or Variable<Matrix>.
        // Because, in these dynamic data type variables, Variable::Zero will create
        // a zero sized vector or matrix, which is useless in the StoreData method.
        // The following method creates correctly sized zero Vector or Matrix
        // accordingly to be assigned to the entities, if they don't have the specified
        // variable in their DataValueContainer.
        const auto& zero = TensorAdaptorUtils::GetZeroValue(*pVariable, this->DataShape());

        ContainerIOUtils::CopyFromContiguousDataArray<data_type>(
            *pContainer, this->ViewData(), r_tensor_shape.data().begin(),
            r_tensor_shape.data().begin() + r_tensor_shape.size(),
            [pVariable, &zero](auto& rEntity) -> auto& {
                if constexpr(DataTypeTraits<data_type>::IsDynamic) {
                    // only dynamic data_type types require the zero
                    // to initialize the uninitialized variables, because
                    // they need to be correctly sized.
                    return rEntity.GetOrCreateValue(*pVariable, zero);
                } else {
                    return rEntity.GetOrCreateValue(*pVariable);
                }
            });
    }, this->mpStorage->GetContainer(), mpVariable);
}

std::string VariableTensorAdaptor::Info() const
{
    std::stringstream info;
    info << "VariableTensorAdaptor:";
    std::visit([&info](auto pVariable) {
        info << " Variable = " << pVariable->Name();
    }, this->mpVariable);
    info << ", " << *(this->mpStorage);
    return info.str();
}

} // namespace Kratos