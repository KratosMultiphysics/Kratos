//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes
#include <cstddef>
#include <type_traits>

namespace Kratos::Internals
{

/**
 * @class StorageView
 * @brief Iterable view over a contiguous array (raw pointer + length).
 * @details Gives a raw storage array the same surface that the uBLAS storage
 * (unbounded_array / std::array) exposes through data(): begin()/end(),
 * size(), operator[] ... plus an implicit conversion to the raw pointer, so
 * both the uBLAS idiom (rV.data().begin()) and the pointer one
 * (Eigen::Map<...>(rV.data(), n), std::memcpy(rV.data(), ...)) compile
 * unchanged. Used by the Eigen-backed dense types for data(), by the CSR
 * wrapper for index1_data()/index2_data()/value_data() and by the amgcl
 * adapters as a range over the CSR arrays.
 * @tparam T The (possibly const-qualified) element type of the viewed array.
 */
template<class T>
class StorageView
{
public:
    using value_type = std::remove_const_t<T>;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using pointer = T*;
    using reference = T&;
    using iterator = T*;
    using const_iterator = const T*;

    /// Constructor wrapping a raw pointer and its length.
    StorageView(T* pData, const std::size_t Size) : mpData(pData), mSize(Size) {}

    /// Iterator to the first element.
    T* begin() const { return mpData; }
    /// Iterator past the last element.
    T* end() const { return mpData + mSize; }
    /// Const iterator to the first element.
    const T* cbegin() const { return mpData; }
    /// Const iterator past the last element.
    const T* cend() const { return mpData + mSize; }
    /// Raw pointer to the first element.
    T* data() const { return mpData; }
    /// Element access by index (any integral index type, so the member is an
    /// exact match and never ambiguous with the built-in pointer subscript
    /// reachable through the implicit pointer conversion).
    template<class TIndexType>
    requires std::is_integral_v<TIndexType>
    T& operator[](const TIndexType Index) const { return mpData[Index]; }
    /// Number of elements viewed.
    std::size_t size() const { return mSize; }
    /// True if the view is empty.
    bool empty() const { return mSize == 0; }
    /// Implicit conversion to the raw storage pointer.
    operator T*() const { return mpData; }

private:
    T* mpData;         /// Pointer to the first element of the viewed storage.
    std::size_t mSize; /// Number of elements viewed.
};

/// Kept for code written against the previous name of the view.
template<class T>
using EigenArrayProxy = StorageView<T>;

} // namespace Kratos::Internals
