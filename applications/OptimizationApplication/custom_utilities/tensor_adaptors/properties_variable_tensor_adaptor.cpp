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
#include <algorithm>
#include <set>
#include <unordered_map>
#include <vector>

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

            const auto p_shape_begin = r_tensor_shape.data().begin();
            const auto p_shape_end = p_shape_begin + r_tensor_shape.size();
            const auto stride = DataTypeTraits<data_type>::Size(p_shape_begin + 1, p_shape_end);
            const auto p_data_begin = this->ViewData().data();

            // group the entities by Properties: entities sharing one must hold equal values, written once
            std::unordered_map<const Properties*, IndexType> first_entity_indices;
            std::vector<IndexType> distinct_entity_indices;
            distinct_entity_indices.reserve(pContainer->size());
            for (IndexType i = 0; i < pContainer->size(); ++i) {
                const auto& r_entity = *(pContainer->begin() + i);
                const auto [p_item, is_new] = first_entity_indices.emplace(&r_entity.GetProperties(), i);
                if (is_new) {
                    distinct_entity_indices.push_back(i);
                } else {
                    const auto p_first_values = p_data_begin + p_item->second * stride;
                    KRATOS_ERROR_IF_NOT(std::equal(p_first_values, p_first_values + stride, p_data_begin + i * stride))
                        << "Entities sharing a properties must have the same " << pVariable->Name()
                        << " value [ properties id = " << r_entity.GetProperties().Id() << ", entity ids = "
                        << (pContainer->begin() + p_item->second)->Id() << " and " << r_entity.Id() << " ].\n";
                }
            }

            // Emplace is not thread safe on one Properties, hence the variable is created serially
            std::vector<data_type*> values(distinct_entity_indices.size());
            if constexpr(DataTypeTraits<data_type>::IsDynamic) {
                // correctly sized zero Vector/Matrix: Variable::Zero() of dynamic types has size zero
                const auto& zero = TensorAdaptorUtils::GetZeroValue(*pVariable, this->DataShape());
                for (IndexType k = 0; k < values.size(); ++k) {
                    values[k] = &(pContainer->begin() + distinct_entity_indices[k])->GetProperties().Emplace(*pVariable, zero);
                }
            } else {
                for (IndexType k = 0; k < values.size(); ++k) {
                    values[k] = &(pContainer->begin() + distinct_entity_indices[k])->GetProperties().Emplace(*pVariable);
                }
            }

            IndexPartition<IndexType>(values.size()).for_each([&](const auto k) {
                DataTypeTraits<data_type>::CopyFromContiguousData(
                    *values[k], p_data_begin + distinct_entity_indices[k] * stride, p_shape_begin + 1, p_shape_end);
            });
        }
    }, this->GetContainer(), mpVariable);
}

std::string PropertiesVariableTensorAdaptor::Info() const
{
    // the base info goes first: it throws for an unloaded (default-constructed) instance, whose mpVariable is null
    const auto base_info = BaseType::Info();
    std::stringstream info;
    info << "PropertiesVariableTensorAdaptor:";
    std::visit([&info](auto pVariable) {
        info << " Variable = " << pVariable->Name();
    }, this->mpVariable);
    info << ", " << base_info;
    return info.str();
}

void PropertiesVariableTensorAdaptor::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
    TensorAdaptorUtils::SaveVariablePointer(rSerializer, "Variable", mpVariable);
}

void PropertiesVariableTensorAdaptor::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
    TensorAdaptorUtils::LoadVariablePointer(rSerializer, "Variable", mpVariable);
}

template KRATOS_API(OPTIMIZATION_APPLICATION) PropertiesVariableTensorAdaptor::PropertiesVariableTensorAdaptor(ModelPart::ConditionsContainerType::Pointer, VariablePointerType);
template KRATOS_API(OPTIMIZATION_APPLICATION) PropertiesVariableTensorAdaptor::PropertiesVariableTensorAdaptor(ModelPart::ConditionsContainerType::Pointer, VariablePointerType, const std::vector<unsigned int>&);
template KRATOS_API(OPTIMIZATION_APPLICATION) PropertiesVariableTensorAdaptor::PropertiesVariableTensorAdaptor(ModelPart::ElementsContainerType::Pointer, VariablePointerType);
template KRATOS_API(OPTIMIZATION_APPLICATION) PropertiesVariableTensorAdaptor::PropertiesVariableTensorAdaptor(ModelPart::ElementsContainerType::Pointer, VariablePointerType, const std::vector<unsigned int>&);

} // namespace Kratos