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
#include <vector>
#include <array> 

// External includes
#include <mdarray.hpp>

// Project includes
#include "containers/mdspan.h" // To reuse aliases like layout_right, extents, etc.

namespace Kratos::Future 
{

/**
 * @brief Provides type aliases and utilities for `mdarray`.
 * @details This header introduces a set of standardized aliases for `mdarray`,
 *          which is an owning multi-dimensional array, contrasting with the non-owning `mdspan`.
 *          It leverages many of the same component aliases as `mdspan.h` (layouts, extents).
 *          Convenient aliases for common `mdarray` types (e.g., 1D and 2D arrays) are defined
 *          to simplify their usage throughout the Kratos codebase.
 */

///@name mdarray Component Aliases (reused from mdspan.h where applicable)
///@{

// Re-expose aliases from mdspan.h for clarity and convenience if not already globally available
// If mdspan.h puts them in Kratos::Future, they are already available.
// For this example, we assume they are in Kratos::Future from mdspan.h
// using Kratos::Future::layout_left;
// using Kratos::Future::layout_right;
// using Kratos::Future::layout_stride;
// template <typename T> using default_accessor = Kratos::Future::default_accessor<T>;
// template <typename TIndexType, std::size_t... TExtents> using extents = Kratos::Future::extents<TIndexType, TExtents...>;
// template <typename TIndexType, std::size_t TRank> using dextents = Kratos::Future::dextents<TIndexType, TRank>;
// inline constexpr auto full_extent = Kratos::Future::full_extent;

///@}
///@name Common mdarray Type Aliases
///@{

/**
 * @brief An owning, multi-dimensional array.
 * @details This is the primary alias for `mdarray`, configured with Kratos defaults.
 *          It owns its data, typically stored in a `std::vector`.
 * @tparam TElementType The type of elements in the array.
 * @tparam TExtents An `extents` object specifying the dimensions.
 * @tparam TLayoutPolicy The memory layout policy (e.g., `layout_right`).
 * @tparam TContainerPolicy The container policy (e.g., `std::vector<TElementType>`).
 */
template <
    typename TElementType,
    typename TExtents,
    typename TLayoutPolicy = Kratos::Future::layout_right, // Default to Kratos standard layout
    typename TContainerPolicy = std::vector<TElementType>
>
using mdarray = MDSPAN_IMPL_STANDARD_NAMESPACE::MDSPAN_IMPL_PROPOSED_NAMESPACE::mdarray<
    TElementType,
    TExtents,
    TLayoutPolicy,
    TContainerPolicy
>;

/**
 * @brief An owning `mdarray` with dynamic extents for a specified rank and a specific container.
 * @tparam T The type of elements in the array.
 * @tparam TRank The number of dimensions (rank).
 * @tparam TLayout The memory layout policy (default is row-major).
 * @tparam TContainer The container type (default is `std::vector<T>`).
 */
template <
    typename T,
    std::size_t TRank,
    typename TLayout = Kratos::Future::layout_right,
    typename TContainer = std::vector<T>
>
using mdarray_dynamic_container = mdarray<T, Kratos::Future::dextents<std::size_t, TRank>, TLayout, TContainer>;

/**
 * @brief An owning `mdarray` with dynamic extents for a specified rank, using `std::vector` as storage.
 * @tparam T The type of elements in the array.
 * @tparam TRank The number of dimensions (rank).
 * @tparam TLayout The memory layout policy (default is row-major).
 */
template <
    typename T,
    std::size_t TRank,
    typename TLayout = Kratos::Future::layout_right
>
using mdarray_dynamic = mdarray_dynamic_container<T, TRank, TLayout, std::vector<T>>;

/**
 * @brief An owning 1D array, conceptually similar to `std::vector` but with mdspan-like interface.
 * @tparam T The type of elements in the array.
 * @tparam TLayout The memory layout policy.
 */
template <typename T, typename TLayout = Kratos::Future::layout_right>
using mdarray_1d = mdarray_dynamic<T, 1, TLayout>;

/**
 * @brief An owning 2D array, conceptually a matrix that owns its data.
 * @tparam T The type of elements in the array.
 * @tparam TLayout The memory layout policy.
 */
template <typename T, typename TLayout = Kratos::Future::layout_right>
using mdarray_2d = mdarray_dynamic<T, 2, TLayout>;

/// @}

} // namespace Kratos::Future
