//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

#pragma once

#include <type_traits>

namespace Kratos {

template <class T>
using BaseType = typename std::remove_cv_t<typename std::remove_reference_t<T>>;

// Primary template that will be used for the fallback
template <class T, typename = void>
struct HasValueType : std::false_type {};

// Specialization for types that have a `value_type` member
template <class T>
struct HasValueType<T, std::void_t<typename T::value_type>> : std::true_type {};

template <typename T>
using IsRValueContainer = std::conjunction<std::is_rvalue_reference<T&&>, HasValueType<T>>;

} // namespace Kratos