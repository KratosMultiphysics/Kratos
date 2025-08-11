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
TensorAdaptor<TDataType>::Storage::Storage(
    ContainerPointerType pContainer,
    const DenseVector<unsigned int>& rShape)
    : mpContainer(pContainer),
      mShape(rShape)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(mShape.empty())
        << "The tensor data shape cannot be empty. It requires at least one dimension representing the number of items in the container [ tensor data = "
        << this->Info() << " ].\n";

    std::visit([this](auto pContainer){
        KRATOS_ERROR_IF_NOT(mShape[0] == pContainer->size())
            << "The value of the first dimension should be equal to the number of items in the pContainer [ container size = "
            << pContainer->size() << ", tensor data = " << this->Info() << " ].\n";
    }, mpContainer);

    // allocate new memory
    mpData = new TDataType[this->Size()];

    KRATOS_CATCH("");
}

template<class TDataType>
TensorAdaptor<TDataType>::Storage::~Storage()
{
    if (mpData) {
        delete[] mpData;
    }
}

template<class TDataType>
typename TensorAdaptor<TDataType>::Storage::Pointer TensorAdaptor<TDataType>::Storage::Copy() const
{
    auto p_tensor_data = Kratos::make_intrusive<Storage>(this->GetContainer(), this->Shape());
    const auto& destination_span = p_tensor_data->ViewData();
    const auto& origin_span = this->ViewData();

    IndexPartition<IndexType>(destination_span.size()).for_each([&destination_span, &origin_span](const auto Index) {
        destination_span[Index] = origin_span[Index];
    });

    return p_tensor_data;
}

template<class TDataType>
Kratos::span<TDataType> TensorAdaptor<TDataType>::Storage::MoveData()
{
    KRATOS_ERROR_IF_NOT(mpData) << "The data is already moved [ " << this->Info() << " ].\n";
    auto p_data = mpData;
    mpData = nullptr;
    return Kratos::span<TDataType>(p_data, p_data + this->Size());
}

template<class TDataType>
Kratos::span<const TDataType>  TensorAdaptor<TDataType>::Storage::ViewData() const
{
    KRATOS_ERROR_IF_NOT(mpData) << "The data is already moved [ " << this->Info() << " ].\n";
    return Kratos::span<const TDataType>(mpData, mpData + this->Size());
}

template<class TDataType>
Kratos::span<TDataType> TensorAdaptor<TDataType>::Storage::ViewData()
{
    KRATOS_ERROR_IF_NOT(mpData) << "The data is already moved [ " << this->Info() << " ].\n";
    return Kratos::span<TDataType>(mpData, mpData + this->Size());
}

template<class TDataType>
DenseVector<unsigned int> TensorAdaptor<TDataType>::Storage::Shape() const
{
    return mShape;
};

template<class TDataType>
DenseVector<unsigned int> TensorAdaptor<TDataType>::Storage::DataShape() const
{
    const auto& shape = this->Shape();
    DenseVector<unsigned int> data_shape(shape.size() - 1);
    std::copy(shape.begin() + 1, shape.end(), data_shape.begin());
    return data_shape;
}

template<class TDataType>
unsigned int TensorAdaptor<TDataType>::Storage::Size() const
{
    return std::accumulate(mShape.data().begin(), mShape.data().end(), 1u, std::multiplies<unsigned int>{});
}

template<class TDataType>
typename TensorAdaptor<TDataType>::Storage::ContainerPointerType TensorAdaptor<TDataType>::Storage::GetContainer() const
{
    return mpContainer;
}

template<class TDataType>
std::string TensorAdaptor<TDataType>::Storage::Info() const
{
    std::stringstream info;
    std::visit([&info, this](auto pContainer) {
        using container_type = std::remove_cv_t<std::decay_t<decltype(*pContainer)>>;
        info << "Storage with " << pContainer->size() << " " << ModelPart::Container<container_type>::GetEntityName() << "(s) with shape = " << this->Shape();
    }, mpContainer);
    return info.str();
}

template<class TDataType>
TensorAdaptor<TDataType>::TensorAdaptor(
    const TensorAdaptor& rOther,
    const bool Copy)
{
    KRATOS_TRY

    if (!Copy) {
        this->mpStorage = rOther.mpStorage;
    } else {
        this->mpStorage = Kratos::make_intrusive<Storage>(rOther.GetContainer(), rOther.Shape());

        const auto& r_current_span = this->ViewData();
        const auto& r_origin_span = rOther.ViewData();
        IndexPartition<IndexType>(rOther.Size()).for_each([&r_current_span, &r_origin_span](const auto Index) {
            r_current_span[Index] = r_origin_span[Index];
        });
    }

    KRATOS_CATCH("");
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
typename TensorAdaptor<TDataType>::Storage::ContainerPointerType TensorAdaptor<TDataType>::GetContainer() const
{
    return mpStorage->GetContainer();
}

template<class TDataType>
Kratos::span<TDataType> TensorAdaptor<TDataType>::MoveData()
{
    return mpStorage->MoveData();
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
    return mpStorage->DataShape();
}

template<class TDataType>
unsigned int TensorAdaptor<TDataType>::Size() const
{
    return mpStorage->Size();
}

template<class TDataType>
std::string TensorAdaptor<TDataType>::Info() const
{
    return mpStorage->Info();
}

// template instantiations
template class TensorAdaptor<bool>;
template class TensorAdaptor<int>;
template class TensorAdaptor<double>;

} // namespace Kratos