//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: HDF5Application/license.txt
//
//  Main author:     Máté Kelemen
//

#pragma once

// Included from "custom_utilities/basic_pipes.h"
#include "custom_utilities/basic_pipes.h" // unnecessary include to get language servers working
#include "includes/exception.h"

// STL includes
#include <cmath>
#include <limits>


namespace Kratos::Pipes {


template <class TVariable>
VariableFromProcessInfo<TVariable>::VariableFromProcessInfo(const Parameters& rParameters)
{
    KRATOS_ERROR_IF_NOT(rParameters.Has("process_info_variable"))
        << "VariableFromProcessInfo expects a Parameters object with \"process_info_variable\", but got:\n"
        << rParameters;
    mVariable = rParameters["process_info_variable"].Get<typename TVariable::Type>();
}


template <class TValue>
IntervalPredicate<TValue>::IntervalPredicate(TValue Begin, TValue End)
    : mInterval(Begin, End)
{
}


template <class TValue>
IntervalPredicate<TValue>::IntervalPredicate(const Parameters& rParameters)
    : mInterval(rParameters)
{
}


template <class TValue>
Modulo<TValue>::Modulo() noexcept
    : mModulo(std::numeric_limits<TValue>::max())
{
}


template <class TValue>
Modulo<TValue>::Modulo(const Parameters& rParameters)
{
    KRATOS_ERROR_IF_NOT(rParameters.Has("mod"))
        << "Modulo expects a Parameters object with a \"mod\" key, but got:\n"
        << rParameters;
    mModulo = rParameters["mod"].Get<TValue>();
}


template <>
inline double Modulo<double>::operator()(double Value) const
{
    return std::fmod(Value, mModulo);
}


} // namespace Kratos::Pipes
