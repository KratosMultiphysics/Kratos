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
#include <unordered_map>
// TODO: Replace with more suited and advanced implementations

// External includes

// Project includes

namespace Kratos {

template<class T1, class T2>
using unordered_map = std::unordered_map<T1, T2>;


} // namespace Kratos