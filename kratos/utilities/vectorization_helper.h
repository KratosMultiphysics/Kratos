//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#pragma once

// System includes
#include <array>
#include <cstddef>
#include <functional>
#include <tuple>
#include <type_traits>
#include <utility>

// External includes

// Project includes
#include "containers/array_1d.h"
#include "containers/nd_data.h"
#include "includes/default_interface.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{

///@addtogroup KratosCore
///@{

namespace detail
{

/**
 * @brief Provides the number of components produced by each call of the vectorized function,
 *        as well as the data type to be stored in the resulting @ref NDData.
 * @details This trait supports scalar results (arithmetic types) as well as fixed size
 *          array like results (std::array, Kratos::array_1d and Kratos::BoundedVector).
 *          Scalar results are stored in an @ref NDData with shape {number_of_entities},
 *          whereas fixed size results are flattened into an @ref NDData with shape
 *          {number_of_entities, components_per_entity}.
 * @note Any other return type is not supported and results in a compilation error.
 */
template<class TDataType, class TEnable = void> struct VectorizationResultInfo;

///@brief Specialization for scalar results.
template<class TDataType>
struct VectorizationResultInfo<TDataType, std::enable_if_t<std::is_arithmetic_v<TDataType>>>
{
    static constexpr unsigned int Size = 1;
    using DataType = TDataType;
};

///@brief Specialization for std::array results.
template<class TDataType, std::size_t TSize>
struct VectorizationResultInfo<std::array<TDataType, TSize>>
{
    static constexpr unsigned int Size = static_cast<unsigned int>(TSize);
    using DataType = TDataType;
};

///@brief Specialization for Kratos::array_1d results.
template<class TDataType, std::size_t TSize>
struct VectorizationResultInfo<array_1d<TDataType, TSize>>
{
    static constexpr unsigned int Size = static_cast<unsigned int>(TSize);
    using DataType = TDataType;
};

///@brief Specialization for Kratos::BoundedVector results.
template<class TDataType, std::size_t TSize>
struct VectorizationResultInfo<BoundedVector<TDataType, TSize>>
{
    static constexpr unsigned int Size = static_cast<unsigned int>(TSize);
    using DataType = TDataType;
};

/**
 * @brief Helper only used in unevaluated contexts to deduce the result type of
 *        invoking @p rFunctor on @p rEntity with the arguments stored in @p rArgs.
 * @details Method calls on the entity (i.e. member function pointers) are
 *          supported by means of @ref std::invoke.
 */
template<class TFunctorType, class TEntityRefType, class TArgsTupleType, std::size_t... TIndices>
auto VectorizationCallHelper(
    TFunctorType& rFunctor,
    TEntityRefType&& rEntity,
    TArgsTupleType& rArgs,
    std::index_sequence<TIndices...>)
    -> decltype(std::invoke(rFunctor, rEntity, std::get<TIndices>(rArgs)...));

///@brief Invokes @p rFunctor on @p rEntity with the arguments stored in @p rArgs.
template<class TFunctorType, class TEntityRefType, class TArgsTupleType, std::size_t... TIndices>
auto InvokeVectorization(
    TFunctorType& rFunctor,
    TEntityRefType&& rEntity,
    TArgsTupleType& rArgs,
    std::index_sequence<TIndices...>)
{
    return std::invoke(rFunctor, rEntity, std::get<TIndices>(rArgs)...);
}

} // namespace detail

/**
 * @brief Applies @p rFunctor to each entity of @p rContainer in parallel, and returns the
 *        collected results in a @ref NDData.
 * @details The function @p rFunctor is invoked once for each of the entries of @p rContainer,
 *          perfectly forwarding @p rArgs to each of the calls. The invocations are distributed
 *          among the threads by means of Kratos parallel utilities (@ref IndexPartition).
 *          The results are stored in a multidimensional @ref NDData:
 *              - scalar results are stored with shape {number_of_entities}
 *              - fixed size results (std::array, array_1d, BoundedVector) are stored with
 *                shape {number_of_entities, components_per_entity}
 *
 *          Note that each of the arguments is used in multiple calls (one per entity).
 *          Hence lvalue arguments are stored as references, whereas arguments passed as
 *          rvalues are copied once and passed to each of the calls by reference.
 *
 * @tparam TContainerType Type of the container to iterate on. Must provide random access iterators.
 * @tparam TFunctorType   Type of the function to be called for each entity.
 * @tparam TArgs          Types of the additional arguments to be forwarded to each call.
 * @param rContainer      Container whose entries @p rFunctor is applied to.
 * @param rFunctor        Function to be called for each entity. Must return either a scalar
 *                        or a fixed size array like type.
 * @param rArgs           Additional arguments perfectly forwarded to each of the calls.
 * @return                Pointer to the @ref NDData with the flattened results.
 */
template<class TContainerType, class TFunctorType, class... TArgs>
[[nodiscard]] auto VectorizationHelper(
    TContainerType& rContainer,
    TFunctorType&& rFunctor,
    TArgs&&... rArgs)
{
    using arg_indices = std::index_sequence_for<TArgs...>;

    // lvalues are stored by reference; rvalues are copied once to be safely used by every call
    using args_storage_type = std::tuple<std::conditional_t<std::is_reference_v<TArgs>, TArgs, std::decay_t<TArgs>>...>;
    args_storage_type stored_args{std::forward<TArgs>(rArgs)...};

    using clean_result_type = std::remove_cv_t<std::remove_reference_t<decltype(detail::VectorizationCallHelper(rFunctor, *rContainer.begin(), stored_args, arg_indices{}))>>;
    using result_info = detail::VectorizationResultInfo<clean_result_type>;
    using nd_data_type = NDData<typename result_info::DataType>;

    DenseVector<unsigned int> shape(result_info::Size == 1u ? 1u : 2u);
    shape[0] = static_cast<unsigned int>(rContainer.size());
    if (result_info::Size != 1u) {
        shape[1] = result_info::Size;
    }

    auto p_nd_data = Kratos::make_shared<nd_data_type>(shape);
    const auto r_data_view = p_nd_data->ViewData();

    const auto it_begin = rContainer.begin();

    IndexPartition<std::size_t>(rContainer.size()).for_each([&](std::size_t Index) {
        if constexpr (result_info::Size == 1u) {
            r_data_view[Index] = detail::InvokeVectorization(
                rFunctor, *(it_begin + Index), stored_args, arg_indices{});
        } else {
            const auto result = detail::InvokeVectorization(
                rFunctor, *(it_begin + Index), stored_args, arg_indices{});
            for (std::size_t i = 0u; i < result_info::Size; ++i) {
                r_data_view[Index * result_info::Size + i] = result[i];
            }
        }
    });

    return p_nd_data;
}

///@}
} // namespace Kratos
