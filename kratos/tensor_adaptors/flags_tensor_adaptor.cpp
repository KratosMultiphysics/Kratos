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
    : mpContainer(pContainer),
      mFlags(rFlags)
{
    this->SetShape(DenseVector<unsigned int>(1, pContainer->size()));
}

FlagsTensorAdaptor::BaseType::Pointer FlagsTensorAdaptor::Clone() const
{
    return std::visit([this](auto pContainer){
        auto p_tensor_adaptor = Kratos::make_intrusive<FlagsTensorAdaptor>(pContainer, this->mFlags);
        IndexPartition<IndexType>(p_tensor_adaptor->Size()).for_each([p_tensor_adaptor, this](const auto Index) {
            p_tensor_adaptor->ViewData()[Index] = this->ViewData()[Index];
        });
        return p_tensor_adaptor;
    }, mpContainer);
}

void FlagsTensorAdaptor::CollectData()
{
    std::visit(
        [this](auto pContainer) {
            const auto& r_tensor_shape = this->Shape();

            using container_type = std::remove_cv_t<std::decay_t<decltype(*pContainer)>>;

            if constexpr (IsInList<container_type,
                                   ModelPart::NodesContainerType,
                                   ModelPart::ConditionsContainerType,
                                   ModelPart::ElementsContainerType>::value) {
                ContainerIOUtils::CopyToContiguousArray<bool>(
                    *pContainer, this->ViewData(), r_tensor_shape.data().begin(),
                    r_tensor_shape.data().begin() + r_tensor_shape.size(),
                    [this](bool& rValue, const auto& rEntity) {
                        rValue = rEntity.Is(this->mFlags);
                    });
            }
            else {
                KRATOS_ERROR
                    << "Flags tensor adaptor can only collect data from nodes, "
                       "conditions or elements containers.";
            }
        },
        mpContainer);
}

void FlagsTensorAdaptor::StoreData()
{
    std::visit(
        [this](auto pContainer) {
            using container_type = std::remove_cv_t<std::decay_t<decltype(*pContainer)>>;

            if constexpr (IsInList<container_type,
                                   ModelPart::NodesContainerType,
                                   ModelPart::ConditionsContainerType,
                                   ModelPart::ElementsContainerType>::value) {
                KRATOS_ERROR_IF_NOT(this->Size() == pContainer->size())
                    << "Size mismatch [ Container size = " << pContainer->size()
                    << ", data span size = " << this->Size() << " ].\n";

                IndexPartition<IndexType>(pContainer->size()).for_each([this, pContainer](const auto Index) {
                    (pContainer->begin() + Index)->Set(this->mFlags, *(this->ViewData().data() + Index));
                });
            }
            else {
                KRATOS_ERROR
                    << "Flags tensor adaptor can only collect data from nodes, "
                       "conditions or elements containers.";
            }
        },
        mpContainer);
}

FlagsTensorAdaptor::ContainerPointerType FlagsTensorAdaptor::GetContainer() const
{
    return this->mpContainer;
}

std::string FlagsTensorAdaptor::Info() const
{
    return std::visit([this](auto pContainer) {
        return TensorAdaptorUtils::Info("FlagsTensorAdaptor ", this->Shape(), *pContainer);
    }, mpContainer);
}

// template instantiations
template KRATOS_API(KRATOS_CORE) FlagsTensorAdaptor::FlagsTensorAdaptor(ModelPart::NodesContainerType::Pointer, const Flags&);
template KRATOS_API(KRATOS_CORE) FlagsTensorAdaptor::FlagsTensorAdaptor(ModelPart::ConditionsContainerType::Pointer, const Flags&);
template KRATOS_API(KRATOS_CORE) FlagsTensorAdaptor::FlagsTensorAdaptor(ModelPart::ElementsContainerType::Pointer, const Flags&);

} // namespace Kratos