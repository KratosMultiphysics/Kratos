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
#include <iostream>
#include <cstddef>

// External includes
#include <mdspan.hpp>

// Project includes

namespace Kratos::Future {

/**
 * @file mdspan.h
 * @brief Provides type aliases and utilities for `std::mdspan`.
 * @details This header introduces a set of standardized aliases for the core
 *          components of `std::mdspan`, such as layouts, accessors, and extents.
 *          It also defines convenient aliases for common `mdspan` types (e.g.,
 *          1D and 2D views) to simplify their usage throughout the Kratos codebase.
 *          A stream insertion operator for easy debugging and inspection of `mdspan`
 *          instances is also included.
 */

///@name mdspan Component Aliases
///@{

/**
 * @brief Alias for a left-to-right (column-major) layout mapping.
 * @see https://en.cppreference.com/w/cpp/container/mdspan/layout_left
 */
using layout_left = MDSPAN_IMPL_STANDARD_NAMESPACE::layout_left;

/**
 * @brief Alias for a right-to-left (row-major) layout mapping.
 * @see https://en.cppreference.com/w/cpp/container/mdspan/layout_right
 */
using layout_right = MDSPAN_IMPL_STANDARD_NAMESPACE::layout_right;

/**
 * @brief Alias for a generalized strided layout mapping.
 * @see https://en.cppreference.com/w/cpp/container/mdspan/layout_stride
 */
using layout_stride = MDSPAN_IMPL_STANDARD_NAMESPACE::layout_stride;

/**
 * @brief Alias for the default accessor policy, which provides element access via `operator[]`.
 * @tparam T The element type being accessed.
 */
template <typename T>
using default_accessor = MDSPAN_IMPL_STANDARD_NAMESPACE::default_accessor<T>;

/**
 * @brief Alias for compile-time extents of an `mdspan`.
 * @tparam TIndexType The type used for indexing.
 * @tparam TExtents A parameter pack of compile-time dimension extents.
 */
template <typename TIndexType, std::size_t... TExtents>
using extents = MDSPAN_IMPL_STANDARD_NAMESPACE::extents<TIndexType, TExtents...>;

/**
 * @brief Alias for run-time (dynamic) extents of an `mdspan`.
 * @tparam TIndexType The type used for indexing.
 * @tparam TRank The number of dimensions (rank) of the `mdspan`.
 */
template <typename TIndexType, std::size_t TRank>
using dextents = MDSPAN_IMPL_STANDARD_NAMESPACE::dextents<TIndexType, TRank>;

/**
 * @brief A constant representing the full extent in an `mdspan`, allowing for dynamic sizing.
 * @details This is used to specify that a dimension should take the full size of the underlying data.
 * @see https://en.cppreference.com/w/cpp/container/mdspan/full_extent
 */
inline constexpr auto full_extent = MDSPAN_IMPL_STANDARD_NAMESPACE::full_extent;

///@}
///@name Common mdspan Type Aliases
///@{

/**
 * @brief A non-owning, multi-dimensional view over a contiguous sequence of objects.
 * @details This is the primary alias for `std::mdspan`, configured with custom
 * layout and accessor policies. It serves as a versatile, non-owning handle
 * to multi-dimensional data.
 * @tparam TElementType The type of elements in the view.
 * @tparam TExtents An `extents` object specifying the dimensions.
 * @tparam TLayoutPolicy The memory layout policy (e.g., `layout_right`).
 * @tparam TAccessorPolicy The policy for accessing elements.
 */
template <
    typename TElementType,
    typename TExtents,
    typename TLayoutPolicy = layout_right,
    typename TAccessorPolicy = default_accessor<TElementType>
>
using mdspan = MDSPAN_IMPL_STANDARD_NAMESPACE::mdspan<
    TElementType,
    TExtents,
    TLayoutPolicy,
    TAccessorPolicy
>;

/**
 * @brief Function template for creating a view of a subset of an mdspan.
 * @details This provides access to the submdspan utility from the implementation namespace.
 * @see https://en.cppreference.com/w/cpp/container/mdspan/submdspan
 */
template <class... TSliceSpecifiers, class TElementType, class TExtents, class TLayout, class TAccessor>
auto submdspan(const mdspan<TElementType, TExtents, TLayout, TAccessor>& source, TSliceSpecifiers... slices)
{
    return MDSPAN_IMPL_STANDARD_NAMESPACE::submdspan(source, slices...);
}

/**
 * @brief A non-owning `mdspan` view with dynamic extents for a specified rank.
 * @tparam T The type of elements in the view.
 * @tparam TRank The number of dimensions (rank).
 * @tparam TLayout The memory layout policy (default is row-major).
 */
template <typename T, std::size_t TRank, typename TLayout = layout_right>
using mdspan_view = mdspan<T, dextents<std::size_t, TRank>, TLayout>;

/**
 * @brief A non-owning 1D view, conceptually equivalent to `std::span`.
 * @tparam T The type of elements in the view.
 */
template <typename T>
using mdspan_1d_view = mdspan_view<T, 1>;

/**
 * @brief A non-owning 2D view, conceptually a view of a matrix.
 * @tparam T The type of elements in the view.
 */
template <typename T>
using mdspan_2d_view = mdspan_view<T, 2>;

/// @}

} // namespace Kratos::Future