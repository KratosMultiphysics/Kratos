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

// System includes
#include "pybind11/pybind11.h"

namespace Kratos::Python {


void AddExpressionIOToPython(pybind11::module& rModule);


} // namespace Kratos::Python
