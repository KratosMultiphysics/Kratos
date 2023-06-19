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

// External includes
#include <pybind11/numpy.h>

// Project includes

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

} // namespace Kratos::Python
