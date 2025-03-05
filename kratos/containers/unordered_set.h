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
// #include <unordered_map> // NOTE: STL implementation

// External includes
#include "ankerl/unordered_dense.h"

// Project includes

namespace Kratos {

template<class T>
using unordered_set = ankerl::unordered_dense::set<T>;
// using unordered_set = std::unordered_set<T>; // NOTE: STL implementation


} // namespace Kratos