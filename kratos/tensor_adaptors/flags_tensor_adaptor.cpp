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
#include <string>

// External includes

// Project includes
#include "containers/flags.h"
#include "tensor_adaptor.h"
#include "tensor_adaptor_utils.h"
#include "utilities/data_type_traits.h"
#include "utilities/container_io_utils.h"

// Include base h
#include "flags_tensor_adaptor.h"

namespace Kratos {

///@name Kratos Classes
///@{

template<class TContainerPointerType>
FlagsTensorAdaptor::FlagsTensorAdaptor(
    TContainerPointerType pContainer,
    const Flags& rFlags)
    : mFlags(rFlags)
{
    this->mpStorage = Kratos::make_intrusive<Storage>(pContainer, DenseVector<unsigned int>(1, pContainer->size()));
}

FlagsTensorAdaptor::FlagsTensorAdaptor(
    const TensorAdaptor& rOther,
    const Flags& rFlags,
    const bool Copy)
    : BaseType(rOther, Copy),
      mFlags(rFlags)
{
    KRATOS_TRY

    const auto& r_tensor_shape = this->mpStorage->Shape();

    KRATOS_ERROR_IF_NOT(r_tensor_shape.size() == 1 && r_tensor_shape[0] == 1)
        << "The data storage within the tensor data is not compatible with the flags "
        << "[ tensor data = " << this->mpStorage->Info() << " ].\n";

    KRATOS_CATCH("");
}

void FlagsTensorAdaptor::Check() const
{
    std::visit([this](auto pContainer) {
        using container_type = BareType<decltype(*pContainer)>;

        if constexpr(IsInList<container_type, ModelPart::NodesContainerType, ModelPart::ConditionsContainerType, ModelPart::ElementsContainerType>) {
            const auto& r_tensor_shape = this->Shape();

            KRATOS_ERROR_IF_NOT(this->Size() == pContainer->size())
                << "Size mismatch [ Container size = " << pContainer->size()
                << ", data span size = " << this->Size() << " ].\n";
        }
    }, this->mpStorage->GetContainer());
}

void FlagsTensorAdaptor::CollectData()
{
    std::visit([this](auto pContainer) {
        using container_type = BareType<decltype(*pContainer)>;

        if constexpr(IsInList<container_type, ModelPart::NodesContainerType, ModelPart::ConditionsContainerType, ModelPart::ElementsContainerType>) {
            const auto& r_tensor_shape = this->Shape();

            KRATOS_ERROR_IF_NOT(this->Size() == pContainer->size())
                << "Size mismatch [ Container size = " << pContainer->size()
                << ", data span size = " << this->Size() << " ].\n";

            ContainerIOUtils::CopyToContiguousArray<bool>(
                *pContainer, this->ViewData(), r_tensor_shape.data().begin(),
                r_tensor_shape.data().begin() + r_tensor_shape.size(),
                [this](bool& rValue, const auto& rEntity) {
                    rValue = rEntity.Is(this->mFlags);
                });
        }
    }, this->mpStorage->GetContainer());
}

void FlagsTensorAdaptor::StoreData()
{
    std::visit([this](auto pContainer) {
        using container_type = BareType<decltype(*pContainer)>;

        if constexpr(IsInList<container_type, ModelPart::NodesContainerType, ModelPart::ConditionsContainerType, ModelPart::ElementsContainerType>) {
            KRATOS_ERROR_IF_NOT(this->Size() == pContainer->size())
                << "Size mismatch [ Container size = " << pContainer->size()
                << ", data span size = " << this->Size() << " ].\n";

            IndexPartition<IndexType>(pContainer->size()).for_each([this, pContainer](const auto Index) {
                (pContainer->begin() + Index)->Set(this->mFlags, *(this->ViewData().data() + Index));
            });
        }
    }, this->mpStorage->GetContainer());
}

std::string FlagsTensorAdaptor::Info() const
{
    std::stringstream info;
    info << "FlagsTensorAdaptor:";
    info << " " << this->mpStorage->Info();
    return info.str();
}

// template instantiations
template KRATOS_API(KRATOS_CORE) FlagsTensorAdaptor::FlagsTensorAdaptor(ModelPart::NodesContainerType::Pointer, const Flags&);
template KRATOS_API(KRATOS_CORE) FlagsTensorAdaptor::FlagsTensorAdaptor(ModelPart::ConditionsContainerType::Pointer, const Flags&);
template KRATOS_API(KRATOS_CORE) FlagsTensorAdaptor::FlagsTensorAdaptor(ModelPart::ElementsContainerType::Pointer, const Flags&);

} // namespace Kratos