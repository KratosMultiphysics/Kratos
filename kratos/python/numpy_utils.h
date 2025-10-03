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

#pragma once

// System includes
#include <numeric>
#include <algorithm>
#include <vector>
#include <iterator>

// External includes
#include <pybind11/numpy.h>

// Project includes
#include "containers/nd_data.h"
#include "utilities/parallel_utilities.h"

namespace Kratos::Python
{

template <class TDataType>
using ContiguousNumpyArray = pybind11::array_t<
    TDataType,
    /*column-major*/ pybind11::array::c_style
>;

template <class TDataType>
pybind11::array_t<TDataType> AllocateNumpyArray(
    const std::size_t NumberOfEntities,
    const std::vector<std::size_t>& rShape)
{
    const std::size_t size = std::accumulate(rShape.begin(),
                                           rShape.end(),
                                           1UL,
                                           [](std::size_t left, std::size_t right) {return left * right;});

    TDataType* array = new TDataType[size * NumberOfEntities];
    pybind11::capsule release(array, [](void* a) {
        delete[] reinterpret_cast<TDataType*>(a);
    });

    // we allocate one additional dimension for the shape to hold the
    // number of entities. So if the shape within kratos is [2, 3] with 70 entities,
    // then the final shape of the numpy array will be [70,2,3] so it can keep the
    // existing shape preserved for each entitiy.
    std::vector<std::size_t> c_shape(rShape.size() + 1);
    c_shape[0] = NumberOfEntities;
    std::copy(rShape.begin(), rShape.end(), c_shape.begin() + 1);

    // we have to allocate number of bytes to be skipped to reach next element
    // in each dimension. So this is calculated backwards.
    std::vector<std::size_t> strides(c_shape.size());
    std::size_t stride_items = 1;
    for (int i = c_shape.size() - 1; i >= 0; --i) {
        strides[i] = sizeof(double) * stride_items;
        stride_items *= c_shape[i];
    }

    return ContiguousNumpyArray<TDataType>(
        c_shape,
        strides,
        array,
        release
    );
}

template <class TDataType>
pybind11::array_t<TDataType> MakeNumpyArray(
    TDataType const* pBegin,
    TDataType const* pEnd,
    const std::vector<std::size_t>& rShape)
{
    auto array = AllocateNumpyArray<TDataType>(rShape);
    KRATOS_ERROR_IF_NOT(std::distance(pBegin, pEnd) == array.size()) << "Size mismatch.";
    std::copy(pBegin,
              pEnd,
              array.mutable_data());
    return array;
}

template<class TDataType, class TPybindArrayType>
bool AssignDataImpl(
    NDData<TDataType>& rDDArray,
    const pybind11::array& rArray)
{
    if (pybind11::isinstance<pybind11::array_t<TPybindArrayType>>(rArray)) {

        KRATOS_ERROR_IF_NOT(rArray.flags() & pybind11::detail::npy_api::constants::NPY_ARRAY_C_CONTIGUOUS_)
            << "Only supports C-style (row-major) arrays from numpy.";

        const auto& cast_array = rArray.cast<pybind11::array_t<TPybindArrayType, pybind11::array::c_style>>();
        auto r_destination_span = rDDArray.ViewData();
        const auto& r_origin_data = cast_array.data();
        IndexPartition<std::size_t>(r_destination_span.size()).for_each([&r_destination_span, &r_origin_data](const auto Index) {
            r_destination_span[Index] = static_cast<TDataType>(r_origin_data[Index]);
        });
        return true;
    } else {
        return false;
    }
}

template<class TDataType, class... TPybindArrayType>
bool AssignData(
    NDData<TDataType>& rDDArray,
    const pybind11::array& rArray)
{
    return (... || AssignDataImpl<TDataType, TPybindArrayType>(rDDArray, rArray));
}

template<class TDataType>
void SetPybindArray(
    NDData<TDataType>& rDDArray,
    const pybind11::array&  rArray)
{
    KRATOS_ERROR_IF(rArray.ndim() == 0)
        << "Passed data is not compatible [ array = "
        << rArray << ", NDData = " << rDDArray << " ].\n";

    std::vector<unsigned int> shape(rArray.ndim());
    std::copy(rArray.shape(), rArray.shape() + rArray.ndim(), shape.begin());

    const auto& r_shape = rDDArray.Shape();

    KRATOS_ERROR_IF_NOT(shape.size() == r_shape.size())
        << "Dimensions mismatch. [ NDData dimensions = " << r_shape.size()
        << ", numpy array dimensions = " << shape.size() << " ].\n";

    for (unsigned int i = 0; i < shape.size(); ++i) {
        KRATOS_ERROR_IF_NOT(r_shape[i] == shape[i])
            << "Shape mismatch. [ NDData shape = " << rDDArray.Shape()
            << ", numpy array shape = " << shape << " ].\n";
    }

    if (!AssignData<
            TDataType,
            bool,
            std::uint8_t,
            std::uint16_t,
            std::uint32_t,
            std::uint64_t,
            std::int8_t,
            std::int16_t,
            std::int32_t,
            std::int64_t,
            float,
            double,
            long double>(rDDArray, rArray))
    {
        KRATOS_ERROR
            << "NDData cannot be assigned an numpy array with \""
            << rArray.dtype() << "\". They can be only set with numpy arrays having following dtypes:"
            << "\n\t numpy.bool"
            << "\n\t numpy.uint8"
            << "\n\t numpy.uint16"
            << "\n\t numpy.uint32"
            << "\n\t numpy.uint64"
            << "\n\t numpy.int8"
            << "\n\t numpy.int16"
            << "\n\t numpy.int32"
            << "\n\t numpy.int64"
            << "\n\t numpy.float32"
            << "\n\t numpy.float64"
            << "\n\t numpy.float128";
    }
}

template<class TDataType>
pybind11::array_t<TDataType, pybind11::array::c_style> GetPybindArray(NDData<TDataType>& rDDArray)
{
    const auto& r_shape = rDDArray.Shape();

    std::vector<std::size_t> c_shape(r_shape.size());
    std::copy(r_shape.begin(), r_shape.end(), c_shape.begin());
    std::vector<std::size_t> strides(c_shape.size());

    std::size_t stride_items = 1;
    for (int i = c_shape.size() - 1; i >= 0; --i) {
        strides[i] = sizeof(TDataType) * stride_items;
        stride_items *= c_shape[i];
    }

    // here we create a raw pointer to the underlying PointerWrapper used in the
    // @p rDDArray. This is done to create an on the fly intrusive_ptr from the
    // raw pointer when creating the capsule, so the PointerWrapper will be kept alive
    // until the capsule releases the memory. This guarantees that, even if the @p rDDArray
    // is killed, the numpy array will be able to function without a problem.
    auto p_raw_data = &*rDDArray.pData();

    pybind11::capsule release(new Kratos::intrusive_ptr<typename NDData<TDataType>::PointerWrapper>(p_raw_data), [](void* a){
        delete static_cast<Kratos::intrusive_ptr<typename NDData<TDataType>::PointerWrapper>*>(a);
    });

    return pybind11::array_t<TDataType, pybind11::array::c_style>(pybind11::buffer_info(
        p_raw_data->Data(),                                 // Pointer to data
        sizeof(TDataType),                                  // Size of one item
        pybind11::format_descriptor<TDataType>::format(),   // Python format descriptor
        c_shape.size(),                                     // Number of dimensions
        c_shape,                                            // Shape of the array
        strides                                             // Strides
    ), release);
}

} // namespace Kratos::Python
