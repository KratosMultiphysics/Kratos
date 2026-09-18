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
#include <set>

// External includes

// Project includes

// Include base h
#include "properties_variable_tensor_adaptor.h"

namespace Kratos {

template<class TContainerPointerType>
PropertiesVariableTensorAdaptor::PropertiesVariableTensorAdaptor(
    TContainerPointerType pContainer,
    VariablePointerType pVariable)
    : mpVariable(pVariable)
{
    this->mpContainer = pContainer;

    std::visit([this, pContainer](auto pVariable) {
        this->mpStorage = Kratos::make_shared<Storage>(
                                TensorAdaptorUtils::GetTensorShape(
                                                *pContainer, *pVariable,
                                                [pVariable](auto& rValue, const auto& rEntity) {
                                                    rValue = rEntity.GetProperties().GetValue(*pVariable);
                                                }));
    }, mpVariable);
}

template<class TContainerPointerType>
PropertiesVariableTensorAdaptor::PropertiesVariableTensorAdaptor(
    TContainerPointerType pContainer,
    VariablePointerType pVariable,
    const std::vector<unsigned int>& rDataShape)
    : mpVariable(pVariable)
{
    this->mpContainer = pContainer;

    std::visit([&rDataShape, this, pContainer](auto pVariable) {
        this->mpStorage = Kratos::make_shared<Storage>(
                                TensorAdaptorUtils::GetTensorShape(
                                *pContainer, *pVariable, rDataShape.data(),
                                rDataShape.data() + rDataShape.size()));
    }, mpVariable);
}

PropertiesVariableTensorAdaptor::PropertiesVariableTensorAdaptor(
    const TensorAdaptor& rOther,
    VariablePointerType pVariable,
    const bool Copy)
    : BaseType(rOther, Copy),
      mpVariable(pVariable)
{
    KRATOS_TRY

    if (!HoldsAlternative<ModelPart::ConditionsContainerType::Pointer,
                          ModelPart::ElementsContainerType::Pointer>::Evaluate(this->GetContainer())) {
        KRATOS_ERROR << "PropertiesVariableTensorAdaptor can only be used with tensor data having condition and element containers "
                     << "[ tensor adaptor = " << rOther << " ].\n";
    }

    // now check whether the given storage is compatible with the variable.
    std::visit([this, &rOther](auto pVariable) {
        using data_type = typename BareType<decltype(*pVariable)>::Type;
        const auto& r_data_shape = this->DataShape();
        KRATOS_ERROR_IF_NOT(DataTypeTraits<data_type>::IsValidShape(r_data_shape.data().begin(), r_data_shape.data().begin() + r_data_shape.size()))
            << "The data storage within the tensor data is not compatible with the " << pVariable->Name()
            << " [ origin data shape = " << rOther.DataShape() << ", tensor adaptor = " << *this << " ].\n";

    }, mpVariable);

    KRATOS_CATCH("");
}

TensorAdaptor<double>::Pointer PropertiesVariableTensorAdaptor::Clone() const
{
    return Kratos::make_shared<PropertiesVariableTensorAdaptor>(*this, mpVariable, true);
}

void PropertiesVariableTensorAdaptor::Check() const
{
    KRATOS_TRY

    BaseType::Check();

    std::visit([](auto pContainer, auto pVariable) {
        using container_type = BareType<decltype(*pContainer)>;

        if constexpr(IsInList<container_type,
                              ModelPart::ConditionsContainerType,
                              ModelPart::ElementsContainerType>) {

            block_for_each(*pContainer, [&pVariable](const auto& rEntity) {
                KRATOS_ERROR_IF_NOT(rEntity.GetProperties().Has(*pVariable))
                    << "The entity with id = " << rEntity.Id() << " does not have the variable " << pVariable->Name() << ".\n";
            });

            // here we have to check for the uniqueness of the variable as well, because
            // otherwise, when StoreData method is called, it can lead to race conditions.
            using variable_type = BareType<decltype(*pVariable)>;

            std::vector<typename variable_type::Type*> memory_locations(pContainer->size());
            IndexPartition(pContainer->size()).for_each([&pVariable, &memory_locations, &pContainer](const auto Index) {
                memory_locations[Index] = &(pContainer->begin() + Index)->GetProperties().GetValue(*pVariable);
            });

            std::set<typename variable_type::Type*> unique_memory_locations(memory_locations.begin(), memory_locations.end());
            KRATOS_ERROR_IF_NOT(unique_memory_locations.size() == memory_locations.size())
                << "The container's entities' properties does not have unique memory locations assigned for variable "
                << pVariable->Name() << " .";
        }
    }, this->GetContainer(), mpVariable);

    KRATOS_CATCH("");
}

void PropertiesVariableTensorAdaptor::CollectData()
{
    std::visit([this](auto pContainer, auto pVariable) {
        using container_type = BareType<decltype(*pContainer)>;

        if constexpr(IsInList<container_type,
                              ModelPart::ConditionsContainerType,
                              ModelPart::ElementsContainerType>) {

            using variable_type = BareType<decltype(*pVariable)>;
            using data_type = typename variable_type::Type;

            const auto& r_tensor_shape = this->Shape();

            KRATOS_ERROR_IF_NOT(r_tensor_shape[0] == pContainer->size())
                << "Underlying container of the tensor data has changed size [ tensor data = "
                << this->Info() << ", container size = " << pContainer->size() << " ].\n";

            ContainerIOUtils::CopyToContiguousArray<data_type>(
                *pContainer, this->ViewData(), r_tensor_shape.data().begin(),
                r_tensor_shape.data().begin() + r_tensor_shape.size(),
                [pVariable](auto& rValue, const auto& rEntity) {
                    rValue = rEntity.GetProperties().GetValue(*pVariable);
                });
        }
    }, this->GetContainer(), mpVariable);
}

void PropertiesVariableTensorAdaptor::StoreData()
{
    std::visit([this](auto pContainer, auto pVariable) {
        using container_type = BareType<decltype(*pContainer)>;

        if constexpr(IsInList<container_type,
                              ModelPart::ConditionsContainerType,
                              ModelPart::ElementsContainerType>) {

            using variable_type = BareType<decltype(*pVariable)>;
            using data_type = typename variable_type::Type;

            const auto& r_tensor_shape = this->Shape();

            KRATOS_ERROR_IF_NOT(r_tensor_shape[0] == pContainer->size())
                << "Underlying container of the tensor data has changed size [ tensor data = "
                << this->Info() << ", container size = " << pContainer->size() << " ].\n";

            if constexpr(DataTypeTraits<data_type>::IsDynamic) {
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
                        return rEntity.GetProperties().Emplace(*pVariable, zero);
                    });
            } else {
                ContainerIOUtils::CopyFromContiguousDataArray<data_type>(
                    *pContainer, this->ViewData(), r_tensor_shape.data().begin(),
                    r_tensor_shape.data().begin() + r_tensor_shape.size(),
                    [pVariable](auto& rEntity) -> auto& {
                        return rEntity.GetProperties().Emplace(*pVariable);
                    });
            }
        }
    }, this->GetContainer(), mpVariable);
}

std::string PropertiesVariableTensorAdaptor::Info() const
{
    std::stringstream info;
    info << "PropertiesVariableTensorAdaptor:";
    std::visit([&info](auto pVariable) {
        info << " Variable = " << pVariable->Name();
    }, this->mpVariable);
    info << ", " << BaseType::Info();
    return info.str();
}

template KRATOS_API(OPTIMIZATION_APPLICATION) PropertiesVariableTensorAdaptor::PropertiesVariableTensorAdaptor(ModelPart::ConditionsContainerType::Pointer, VariablePointerType);
template KRATOS_API(OPTIMIZATION_APPLICATION) PropertiesVariableTensorAdaptor::PropertiesVariableTensorAdaptor(ModelPart::ConditionsContainerType::Pointer, VariablePointerType, const std::vector<unsigned int>&);
template KRATOS_API(OPTIMIZATION_APPLICATION) PropertiesVariableTensorAdaptor::PropertiesVariableTensorAdaptor(ModelPart::ElementsContainerType::Pointer, VariablePointerType);
template KRATOS_API(OPTIMIZATION_APPLICATION) PropertiesVariableTensorAdaptor::PropertiesVariableTensorAdaptor(ModelPart::ElementsContainerType::Pointer, VariablePointerType, const std::vector<unsigned int>&);

} // namespace Kratos