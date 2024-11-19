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

template<class T1, class T2>
using unordered_map = ankerl::unordered_dense::map<T1, T2>;
// using unordered_map = std::unordered_map<T1, T2>; // NOTE: STL implementation


} // namespace Kratos