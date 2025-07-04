//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A Roig

#pragma once

// System includes
#include <unordered_set>

// External includes

// Project includes
#include "solving_strategies/builder_and_solvers/custom_allocators/intel_tbb.h"

namespace Kratos {

// This namespace containers wrappers to std containers globally modifying its allocators
// to avoid having to include all possible allocators in every file where a custom allocator
// is needed.

// Usually we would preffer to have a KratosContainers class but namespace is needed
// unless we want to include all reference container in the same file, which is not correct
namespace Containers {

    // Unordered Set
    template<
        class TDataType,
        class THashType  = std::hash<TDataType>,
        class TEqualType = std::equal_to<TDataType>,
#ifdef KRATOS_INTEL_TBB
        class TAllocator = tbb::scalable_allocator<TDataType>
#else
        class TAllocator = std::allocator<TDataType>
#endif
    >
    using unordered_set = std::unordered_set<TDataType, THashType, TEqualType, TAllocator>;
}
}