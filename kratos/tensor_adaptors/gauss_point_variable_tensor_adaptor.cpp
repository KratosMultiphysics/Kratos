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

template<class TContainerPointerType>
GaussPointVariableTensorAdaptor::GaussPointVariableTensorAdaptor(
    TContainerPointerType pContainer,
    VariablePointerType pVariable,
    ProcessInfo::Pointer pProcessInfo)
    : mpVariable(pVariable),
      mpProcessInfo(pProcessInfo)
{
    using container_type = BareType<decltype(*pContainer)>;
    using entity_type = typename container_type::value_type;

    std::visit([this, pContainer](auto pVariable) {
        using variable_type = BareType<decltype(*pVariable)>;
        using data_type = typename variable_type::Type;

        this->mpStorage = Kratos::make_intrusive<TensorData<double>>(
            pContainer,
            TensorAdaptorUtils::GetTensorShape<std::vector<data_type>>(
                *pContainer, [pVariable, this](auto& rValues, const auto& rEntity) {
                    const_cast<entity_type&>(rEntity).CalculateOnIntegrationPoints(*pVariable, rValues, *(this->mpProcessInfo));
                }));
    }, mpVariable);
}

void GaussPointVariableTensorAdaptor::CollectData()
{
    std::visit([this](auto pContainer, auto pVariable) {
        using container_type = BareType<decltype(*pContainer)>;

        if constexpr(IsInList<container_type, ModelPart::ConditionsContainerType, ModelPart::ElementsContainerType>) {
            using entity_type = typename container_type::value_type;
            using variable_type = BareType<decltype(*pVariable)>;
            using data_type = typename variable_type::Type;

            const auto& tensor_shape = this->Shape();

            ContainerIOUtils::CopyToContiguousArray<std::vector<data_type>>(
                *pContainer, this->ViewData(), tensor_shape.data().begin(),
                tensor_shape.data().begin() + tensor_shape.size(),
                [pVariable, this](auto& rValues, const auto& rEntity) {
                    const_cast<entity_type&>(rEntity).CalculateOnIntegrationPoints(
                        *pVariable, rValues, *(this->mpProcessInfo));
                });
        }
    }, this->mpStorage->GetContainer(), mpVariable);
}

void GaussPointVariableTensorAdaptor::StoreData()
{
    KRATOS_ERROR << "Storing gauss point data is not supported.";
}

std::string GaussPointVariableTensorAdaptor::Info() const
{
    std::stringstream info;
    info << "GaussPointVariableTensorAdaptor:";
    std::visit([&info](auto pVariable) {
        info << " Variable = " << pVariable->Name();
    }, this->mpVariable);
    info << ", " << *(this->mpStorage);
    return info.str();
}

// template instantiations
template KRATOS_API(KRATOS_CORE) GaussPointVariableTensorAdaptor::GaussPointVariableTensorAdaptor(ModelPart::ConditionsContainerType::Pointer, VariablePointerType, ProcessInfo::Pointer);
template KRATOS_API(KRATOS_CORE) GaussPointVariableTensorAdaptor::GaussPointVariableTensorAdaptor(ModelPart::ElementsContainerType::Pointer, VariablePointerType, ProcessInfo::Pointer);

} // namespace Kratos