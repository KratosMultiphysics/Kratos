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

// --- Utility Includes ---
#include "utilities/profiler.h" // <== help the language server

// --- STL Includes ---
#include <string>
#include <chrono>


namespace Kratos::Internals {


template <class T>
Profiler<T>::Scope::Scope(Profiler::Item& rItem)
    : mrItem(rItem),
      mBegin(Clock::now())
{
    ++mrItem.mCallCount;
    ++mrItem.mRecursionLevel;
}


template <class T>
Profiler<T>::Scope::~Scope()
{
    if (!--mrItem.mRecursionLevel)
        mrItem.mTime += std::chrono::duration_cast<Profiler::TimeUnit>(Clock::now() - mBegin);
}


template <class T>
typename Profiler<T>::Scope Profiler<T>::Profile(Item& rItem)
{
    return Scope(rItem);
}


} // namespace Kratos::Internals
