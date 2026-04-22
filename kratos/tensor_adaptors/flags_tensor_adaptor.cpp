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
    : BaseType(pContainer,  Kratos::make_shared<Storage>(DenseVector<unsigned int>(1, pContainer->size()))),
      mFlags(rFlags)
{
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

    if (!HoldsAlternative<ModelPart::NodesContainerType::Pointer,
                          ModelPart::ConditionsContainerType::Pointer,
                          ModelPart::ElementsContainerType::Pointer>::Evaluate(this->GetContainer())) {
        KRATOS_ERROR << "FlagsTensorAdaptor can only be used with tensor data having nodal, condition or element containers "
                     << "[ tensor adaptor = " << rOther << " ].\n";
    }

    std::visit([&r_tensor_shape, &rOther](auto pContainer){
        KRATOS_ERROR_IF_NOT(r_tensor_shape.size() == 1 && r_tensor_shape[0] == pContainer->size())
            << "The data storage within the tensor data is not compatible with the flags "
            << "[ tensor adaptor = " << rOther << " ].\n";
    }, this->GetContainer());

    KRATOS_CATCH("");
}

void FlagsTensorAdaptor::Check() const
{
    KRATOS_TRY

    std::visit([this](auto pContainer) {
        using container_type = BareType<decltype(*pContainer)>;

        if constexpr(IsInList<container_type, ModelPart::NodesContainerType, ModelPart::ConditionsContainerType, ModelPart::ElementsContainerType>) {
            KRATOS_ERROR_IF_NOT(this->Size() == pContainer->size())
                << "Size mismatch [ Container size = " << pContainer->size()
                << ", data span size = " << this->Size() << " ].\n";
        }
    }, this->GetContainer());

    KRATOS_CATCH("");
}

TensorAdaptor<int>::Pointer FlagsTensorAdaptor::Clone() const
{
    return Kratos::make_shared<FlagsTensorAdaptor>(*this);
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

            ContainerIOUtils::CopyToContiguousArray<int>(
                *pContainer, this->ViewData(), r_tensor_shape.data().begin(),
                r_tensor_shape.data().begin() + r_tensor_shape.size(),
                [this](int& rValue, const auto& rEntity) {
                    rValue = (rEntity.IsDefined(this->mFlags) ? rEntity.Is(this->mFlags) : -1);
                });
        }
    }, this->GetContainer());
}

void FlagsTensorAdaptor::StoreData()
{
    std::visit([this](auto pContainer) {
        using container_type = BareType<decltype(*pContainer)>;

        if constexpr(IsInList<container_type, ModelPart::NodesContainerType, ModelPart::ConditionsContainerType, ModelPart::ElementsContainerType>) {
            KRATOS_ERROR_IF_NOT(this->Size() == pContainer->size())
                << "Size mismatch [ Container size = " << pContainer->size()
                << ", data span size = " << this->Size() << " ].\n";

            auto span_itr = this->ViewData().data();

            IndexPartition<IndexType>(pContainer->size()).for_each([this, span_itr, pContainer](const auto Index) {
                const int value = *(span_itr + Index);
                auto& r_entity = *(pContainer->begin() + Index);
                switch (value) {
                    case 0:
                    case 1:
                        r_entity.Set(this->mFlags, value);
                        break;
                    case -1:
                        r_entity.Reset(this->mFlags);
                        break;
                    default:
                        KRATOS_ERROR << "Invalid flag status = " << value << " for " << ModelPart::Container<container_type>::GetEntityName()
                                     << " with id " <<  r_entity.Id() << ". Flag tensor adaptor's internal storage should only have following values:"
                                     << "\n\t -1 - The corresponding flag is not defined."
                                     << "\n\t  0 - The corresponding flag is defined and false."
                                     << "\n\t  1 - The corresponding flag is defined and true.";
                        break;
                }
            });
        }
    }, this->GetContainer());
}

std::string FlagsTensorAdaptor::Info() const
{
    std::stringstream info;
    info << "FlagsTensorAdaptor:";
    info << " " << BaseType::Info();
    return info.str();
}

// template instantiations
template KRATOS_API(KRATOS_CORE) FlagsTensorAdaptor::FlagsTensorAdaptor(ModelPart::NodesContainerType::Pointer, const Flags&);
template KRATOS_API(KRATOS_CORE) FlagsTensorAdaptor::FlagsTensorAdaptor(ModelPart::ConditionsContainerType::Pointer, const Flags&);
template KRATOS_API(KRATOS_CORE) FlagsTensorAdaptor::FlagsTensorAdaptor(ModelPart::ElementsContainerType::Pointer, const Flags&);

} // namespace Kratos