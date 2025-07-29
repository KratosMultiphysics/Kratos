//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

#pragma once

// Project includes
#include "utilities/comparison_impl.hpp" // IntegerComparison, FloatComparison
#include "includes/kratos_export_api.h" // KRATOS_API


namespace Kratos {


/// @brief Unspecialized base class for comparison operators to be specialized for specific types.
template <class T>
struct KRATOS_API(KRATOS_CORE) Comparison {
    static_assert(std::is_same_v<T,void>, "attempting to instantiate Comparison for unsupported type");
    struct Equal {};
    struct Less {};
};


template <>
struct KRATOS_API(KRATOS_CORE) Comparison<int> : public Impl::IntegerComparison<int> {
    using Impl::IntegerComparison<int>::IntegerComparison;
}; // Comparison<int>


template <>
struct KRATOS_API(KRATOS_CORE) Comparison<unsigned> : public Impl::IntegerComparison<unsigned> {
    using Impl::IntegerComparison<unsigned>::IntegerComparison;
}; // Comparison<unsigned>


template <>
struct KRATOS_API(KRATOS_CORE) Comparison<std::size_t> : public Impl::IntegerComparison<std::size_t> {
    using Impl::IntegerComparison<std::size_t>::IntegerComparison;
}; // Comparison<std::size_t>


/// @copydoc Impl::FloatComparison
template <>
struct KRATOS_API(KRATOS_CORE) Comparison<float> : public Impl::FloatComparison<float> {
    using Impl::FloatComparison<float>::FloatComparison;
}; // Comparison<float>


/// @copydoc Impl::FloatComparison
template <>
struct KRATOS_API(KRATOS_CORE) Comparison<double> : public Impl::FloatComparison<double> {
    using Impl::FloatComparison<double>::FloatComparison;
}; // Comparison<double>


} // namespace Kratos
