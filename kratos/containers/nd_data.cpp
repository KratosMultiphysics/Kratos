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
#include <numeric>
#include <sstream>

// External	includes

// Project includes
#include "includes/define.h"
#include "utilities/parallel_utilities.h"

// Include base h
#include "nd_data.h"

namespace Kratos {

template<class TDataType>
NDData<TDataType>::NDData(const DenseVector<unsigned int>& rShape)
    : mShape(rShape)
{
    // allocate new memory
    mpData = Kratos::make_intrusive<PointerWrapper>(new TDataType[this->Size()], true);
}

template<class TDataType>
NDData<TDataType>::NDData(
    const DenseVector<unsigned int>& rShape,
    const TDataType Value)
    : NDData(rShape)
{
    auto span = this->ViewData();
    IndexPartition<IndexType>(this->Size()).for_each([&span, Value](const auto Index) {
        span[Index] = Value;
    });
}

template<class TDataType>
NDData<TDataType>::NDData(
    TDataType * pData,
    const DenseVector<unsigned int>& rShape,
    const bool Copy)
    : mShape(rShape)
{
    if (!Copy) {
        mpData = Kratos::make_intrusive<PointerWrapper>(pData, false);
    } else {
        mpData = Kratos::make_intrusive<PointerWrapper>(new TDataType[this->Size()], true);

        const auto span = this->ViewData();
        IndexPartition<IndexType>(this->Size()).for_each([&span, pData](const auto Index) {
            span[Index] = pData[Index];
        });
    }
}

template<class TDataType>
NDData<TDataType>::NDData(
    TDataType const * pData,
    const DenseVector<unsigned int>& rShape)
    : mShape(rShape)
{
    mpData = Kratos::make_intrusive<PointerWrapper>(new TDataType[this->Size()], true);

    const auto span = this->ViewData();
    IndexPartition<IndexType>(this->Size()).for_each([&span, pData](const auto Index) {
        span[Index] = pData[Index];
    });
}

template<class TDataType>
NDData<TDataType>::NDData(const NDData& rOther)
    : NDData(rOther.ViewData().data(), rOther.Shape())
{
}

template<class TDataType>
Kratos::span<const TDataType> NDData<TDataType>::ViewData() const
{
    return Kratos::span<const TDataType>(this->mpData->Data(), this->mpData->Data() + this->Size());
}

template<class TDataType>
Kratos::span<TDataType> NDData<TDataType>::ViewData()
{
    return Kratos::span<TDataType>(this->mpData->Data(), this->mpData->Data() + this->Size());
}

template<class TDataType>
typename NDData<TDataType>::PointerWrapper::Pointer NDData<TDataType>::pData() const
{
    return mpData;
}

template<class TDataType>
DenseVector<unsigned int> NDData<TDataType>::Shape() const
{
    return mShape;
}

template<class TDataType>
unsigned int NDData<TDataType>::Size() const
{
    return std::accumulate(mShape.data().begin(), mShape.data().end(), 1u, std::multiplies<unsigned int>{});
}

template<class TDataType>
std::string NDData<TDataType>::Info() const
{
    std::stringstream info;
    info << "NDData with shape = " << this->Shape() << ".";
    return info.str();
}

// template instantiations
template class NDData<unsigned char>; // We have to use the unsigned char, because numpy does not have proper bindings for std::uint8_t.
template class NDData<bool>;
template class NDData<int>;
template class NDData<double>;

} // namespace Kratos