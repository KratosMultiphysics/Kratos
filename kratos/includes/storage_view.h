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
 * @brief Iterable view over a contiguous array (raw pointer and length)
 * @details Gives a raw storage array the same surface that the uBLAS storage (unbounded_array / std::array) exposes
 * through data(): begin()/end(), size(), operator[] ... plus an implicit conversion to the raw pointer, so both the
 * uBLAS idiom (rV.data().begin()) and the pointer one (Eigen::Map<...>(rV.data(), n), std::memcpy(rV.data(), ...))
 * compile unchanged. Used by the Eigen-backed dense types for data(), by the CSR wrapper for
 * index1_data()/index2_data()/value_data() and by the amgcl adapters as a range over the CSR arrays.
 * @tparam T The (possibly const-qualified) element type of the viewed array.
 */
template<class T>
class StorageView
{
public:
    /// Element type of the view. The viewed element type without const qualification.
    using value_type = std::remove_const_t<T>;
    /// Type used for element counts. std::size_t, as in the uBLAS containers.
    using size_type = std::size_t;
    /// Signed type of the distance between two iterators. std::ptrdiff_t.
    using difference_type = std::ptrdiff_t;
    /// Pointer to an element of the view. T*, the type returned by data().
    using pointer = T*;
    /// Reference to an element of the view. T&, the type returned by operator[].
    using reference = T&;
    /// Iterator over the viewed elements. T*, the iterator type of begin() and end().
    using iterator = T*;
    /// Const iterator over the viewed elements. const T*, the iterator type of cbegin() and cend().
    using const_iterator = const T*;

    /**
     * @brief Constructs the view from a raw pointer and a length
     * @details Stores the pointer and the length without copying the storage, so the viewed storage must outlive the
     * view.
     * @param pData Pointer to the first element of the viewed storage.
     * @param Size Number of elements viewed.
     */
    StorageView(T* pData, const std::size_t Size) : mpData(pData), mSize(Size) {}

    /**
     * @brief Iterator to the first element
     * @details Returns the raw pointer to the first element, so the view can be used in range-based loops.
     * @return Pointer to the first element.
     */
    T* begin() const { return mpData; }
    /**
     * @brief Iterator past the last element
     * @details Points to mpData + size(), so that [begin(), end()) covers the viewed elements.
     * @return Pointer past the last element.
     */
    T* end() const { return mpData + mSize; }
    /**
     * @brief Const iterator to the first element
     * @details Same as begin(), with a const pointer.
     * @return Const pointer to the first element.
     */
    const T* cbegin() const { return mpData; }
    /**
     * @brief Const iterator past the last element
     * @details Same as end(), with a const pointer.
     * @return Const pointer past the last element.
     */
    const T* cend() const { return mpData + mSize; }
    /**
     * @brief Raw pointer to the first element
     * @details Returns the pointer to the viewed storage, as the uBLAS data() accessor does.
     * @return Pointer to the first element of the viewed storage.
     */
    T* data() const { return mpData; }
    /**
     * @brief Element access by index
     * @details Accepts any integral index type. Being a template, the member is an exact match, so it is never
     * ambiguous with the built-in pointer subscript reachable through the implicit pointer conversion.
     * @tparam TIndexType Integral type of the index.
     * @param Index Position of the element.
     * @return Reference to the element at position Index.
     */
    template<class TIndexType>
    requires std::is_integral_v<TIndexType>
    T& operator[](const TIndexType Index) const { return mpData[Index]; }
    /**
     * @brief Number of elements viewed
     * @details Returns the length of the view.
     * @return Number of elements viewed.
     */
    std::size_t size() const { return mSize; }
    /**
     * @brief True if the view is empty
     * @details Returns true when the view holds no elements.
     * @return True if the view holds no elements.
     */
    bool empty() const { return mSize == 0; }
    /**
     * @brief Implicit conversion to the raw storage pointer
     * @details Returns the pointer to the first element, so the view can be passed where a raw pointer is expected.
     * @return Pointer to the first element, as the raw storage pointer.
     */
    operator T*() const { return mpData; }

private:
    /**
     * @brief Pointer to the first element of the viewed storage
     * @details Not owned by the view.
     */
    T* mpData;
    /**
     * @brief Number of elements viewed
     * @details Stored as a std::size_t.
     */
    std::size_t mSize;
};

/// Previous name of StorageView. Alias template kept for code written against the previous name of the view. It is
/// exactly StorageView<T>.
template<class T>
using EigenArrayProxy = StorageView<T>;

} // namespace Kratos::Internals
