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

        this->mpStorage = Kratos::make_intrusive<Storage>(
            pContainer,
            TensorAdaptorUtils::GetTensorShape<std::vector<data_type>>(
                *pContainer, [pVariable, this](auto& rValues, const auto& rEntity) {
                    const_cast<entity_type&>(rEntity).CalculateOnIntegrationPoints(*pVariable, rValues, *(this->mpProcessInfo));
                }));
    }, mpVariable);
}

GaussPointVariableTensorAdaptor::GaussPointVariableTensorAdaptor(
    const TensorAdaptor& rOther,
    VariablePointerType pVariable,
    ProcessInfo::Pointer pProcessInfo,
    const bool Copy)
    : BaseType(rOther, Copy),
      mpVariable(pVariable),
      mpProcessInfo(pProcessInfo)
{
    KRATOS_TRY

    // now check whether the given storage is compatible with the variable.
    std::visit([this, &rOther](auto pVariable) {
        using data_type = typename BareType<decltype(*pVariable)>::Type;
        const auto& r_data_shape = this->mpStorage->DataShape();
        KRATOS_ERROR_IF_NOT(DataTypeTraits<std::vector<data_type>>::IsValidShape(r_data_shape.data().begin(), r_data_shape.data().begin() + r_data_shape.size()))
            << "The data storage within the tensor data is not compatible with the " << pVariable->Name()
            << " [ origin data shape = " << rOther.DataShape() << ", tensor adaptor = " << *this << " ].\n";

    }, mpVariable);

    KRATOS_CATCH("");
}

GaussPointVariableTensorAdaptor::GaussPointVariableTensorAdaptor(
    const GaussPointVariableTensorAdaptor& rOther,
    const bool Copy)
    : GaussPointVariableTensorAdaptor(rOther, rOther.mpVariable, rOther.mpProcessInfo, Copy)
{
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

            KRATOS_ERROR_IF_NOT(tensor_shape[0] == pContainer->size())
                << "Underlying container of the tensor data has changed size [ tensor data = "
                << this->mpStorage->Info() << ", container size = " << pContainer->size() << " ].\n";

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
    info << ", " << this->mpStorage->Info();
    return info.str();
}

// template instantiations
template KRATOS_API(KRATOS_CORE) GaussPointVariableTensorAdaptor::GaussPointVariableTensorAdaptor(ModelPart::ConditionsContainerType::Pointer, VariablePointerType, ProcessInfo::Pointer);
template KRATOS_API(KRATOS_CORE) GaussPointVariableTensorAdaptor::GaussPointVariableTensorAdaptor(ModelPart::ElementsContainerType::Pointer, VariablePointerType, ProcessInfo::Pointer);

} // namespace Kratos