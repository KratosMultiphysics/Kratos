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
#include "utilities/profiler.h" // <== help the language server

// System includes
#include <chrono>


namespace Kratos::Internals {


template <class T>
Profiler<T>::Scope::Scope(Profiler::Item& rItem)
    : Scope(rItem, Clock::now())
{
}


template <class T>
Profiler<T>::Scope::Scope(Profiler::Item& rItem, std::chrono::high_resolution_clock::time_point Begin)
    : mrItem(rItem),
      mBegin(Begin)
{
    ++mrItem.mCallCount;
    ++mrItem.mRecursionLevel;
}


template <class T>
Profiler<T>::Scope::~Scope()
{
    if (!--mrItem.mRecursionLevel) {
        const auto duration = std::chrono::duration_cast<Profiler::TimeUnit>(Clock::now() - mBegin);
        mrItem.mCumulative += duration;
        mrItem.mMin = std::min(mrItem.mMin, duration);
        mrItem.mMax = std::max(mrItem.mMax, duration);
    }
}


template <class T>
typename Profiler<T>::Scope Profiler<T>::Profile(Item& rItem)
{
    return Scope(rItem);
}


} // namespace Kratos::Internals
