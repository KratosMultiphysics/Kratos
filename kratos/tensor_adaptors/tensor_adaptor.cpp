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
#include <atomic>
#include <variant>
#include <numeric>

// External includes
#include <span/span.hpp>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "intrusive_ptr/intrusive_ptr.hpp"
#include "utilities/parallel_utilities.h"

// Include base h
#include "tensor_adaptor.h"

namespace Kratos {

template<class TDataType>
TensorAdaptor<TDataType>::TensorAdaptor(
    ContainerPointerType pContainer,
    typename NDData<TDataType>::Pointer pData,
    const bool Copy)
{
    if (!Copy) {
        this->mpStorage = pData;
    } else {
        this->mpStorage = Kratos::make_shared<NDData<TDataType>>(pData->ViewData().data(), pData->Shape(), Copy);
    }
    this->mpContainer = pContainer;
}

template<class TDataType>
TensorAdaptor<TDataType>::TensorAdaptor(
    const TensorAdaptor& rOther,
    const bool Copy)
{
    if (!Copy) {
        this->mpStorage = rOther.mpStorage;
    } else {
        this->mpStorage = Kratos::make_shared<NDData<TDataType>>(rOther.ViewData().data(), rOther.Shape());
    }
    this->mpContainer = rOther.mpContainer;
}

template<class TDataType>
TensorAdaptor<TDataType>::Pointer TensorAdaptor<TDataType>::Clone() const
{
    return Kratos::make_shared<TensorAdaptor<TDataType>>(*this);
}

template<class TDataType>
void TensorAdaptor<TDataType>::Check() const
{
    KRATOS_TRY

    std::visit([this](auto pContainer) {
        const auto& r_tensor_shape = this->Shape();
        KRATOS_ERROR_IF_NOT(r_tensor_shape[0] == pContainer->size())
            << "The underlying container size mismatch with the first dimension of the tensor shape "
            << "[ tensor = " << *this << " ].\n";
    }, this->GetContainer());

    KRATOS_CATCH("");
}

template<class TDataType>
void TensorAdaptor<TDataType>::CollectData()
{
    KRATOS_ERROR << "Calling TensorAdaptor::CollectData method. This base class can only be used for data storage, not to collect or store data.";
}

template<class TDataType>
void TensorAdaptor<TDataType>::StoreData()
{
    KRATOS_ERROR << "Calling TensorAdaptor::StoreData method. This base class can only be used for data storage, not to collect or store data.";
}

template<class TDataType>
typename TensorAdaptor<TDataType>::ContainerPointerType TensorAdaptor<TDataType>::GetContainer() const
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(mpContainer.has_value())
        << "Tensor adaptor is created without a container [ tensor adaptor = " << this->Info() << " ].\n";
    return mpContainer.value();

    KRATOS_CATCH("");
}

template<class TDataType>
bool TensorAdaptor<TDataType>::HasContainer() const
{
    KRATOS_TRY

    return mpContainer.has_value();

    KRATOS_CATCH("");
}

template<class TDataType>
Kratos::span<const TDataType> TensorAdaptor<TDataType>::ViewData() const
{
    const auto& storage = *mpStorage;
    return storage.ViewData();
}

template<class TDataType>
Kratos::span<TDataType> TensorAdaptor<TDataType>::ViewData()
{
    return mpStorage->ViewData();
}

template<class TDataType>
DenseVector<unsigned int> TensorAdaptor<TDataType>::Shape() const
{
    return mpStorage->Shape();
}

template<class TDataType>
DenseVector<unsigned int> TensorAdaptor<TDataType>::DataShape() const
{
    const auto& r_shape = this->Shape();
    DenseVector<unsigned int> data_shape(r_shape.size() - 1);
    std::copy(r_shape.begin() + 1, r_shape.end(), data_shape.begin());
    return data_shape;
}

template<class TDataType>
unsigned int TensorAdaptor<TDataType>::Size() const
{
    return mpStorage->Size();
}

template<class TDataType>
typename TensorAdaptor<TDataType>::Storage::Pointer TensorAdaptor<TDataType>::pGetStorage()
{
    return mpStorage;
}

template<class TDataType>
std::string TensorAdaptor<TDataType>::Info() const
{
    std::stringstream info;

    if (mpContainer.has_value()) {
        std::visit([&info, this](auto pContainer) {
            using container_type = std::remove_cv_t<std::decay_t<decltype(*pContainer)>>;
            info << "Storage with " << pContainer->size() << " " << ModelPart::Container<container_type>::GetEntityName() << "(s) with shape = " << this->Shape();
        }, mpContainer.value());
    } else {
        info << "Storage with with shape = " << this->Shape();
    }
    return info.str();
}

// template instantiations
template class TensorAdaptor<bool>;
template class TensorAdaptor<int>;
template class TensorAdaptor<double>;

} // namespace Kratos