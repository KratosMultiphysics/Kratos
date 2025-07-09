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
#include "tensor_adaptor.h"
#include "tensor_adaptor_utils.h"
#include "utilities/container_io_utils.h"

// Include base h
#include "flags_tensor_adaptor.h"

namespace Kratos {

template<class TContainerPointerType>
FlagsTensorAdaptor::FlagsTensorAdaptor(
    TContainerPointerType pContainer,
    const Flags& rFlags)
    : mpIO(Kratos::make_shared<FlagsIO>(rFlags)),
      mpContainer(pContainer)
{
    // set the shape
    this->mShape.resize(1);
    this->mShape[0] = pContainer->size();

    TensorAdaptorUtils::InitializeData(*pContainer, *mpIO, this->Shape(), this->mData);
}

void FlagsTensorAdaptor::CollectData()
{
    std::visit([this](auto pContainer) {
        TensorAdaptorUtils::CollectData(*pContainer, *this->mpIO, this->Shape(), this->mData);
    }, this->mpContainer);
}

void FlagsTensorAdaptor::StoreData()
{
    std::visit([this](auto pContainer) {
        TensorAdaptorUtils::StoreData(*pContainer, *this->mpIO, this->Shape(), this->mData);
    }, this->mpContainer);
}

FlagsTensorAdaptor::ContainerType FlagsTensorAdaptor::GetContainer() const
{
    return mpContainer;
}

std::string FlagsTensorAdaptor::Info() const
{
    return std::visit([this](auto pContainer) {
        return TensorAdaptorUtils::Info(*pContainer, *mpIO, this->Shape());
    }, mpContainer);
}

// template instantiations
template FlagsTensorAdaptor::FlagsTensorAdaptor(ModelPart::NodesContainerType::Pointer, const Flags&);
template FlagsTensorAdaptor::FlagsTensorAdaptor(ModelPart::ConditionsContainerType::Pointer, const Flags&);
template FlagsTensorAdaptor::FlagsTensorAdaptor(ModelPart::ElementsContainerType::Pointer, const Flags&);

} // namespace Kratos