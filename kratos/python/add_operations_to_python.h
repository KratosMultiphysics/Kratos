//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos Roig
//                   Ruben Zorrilla
//

#pragma once

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes

namespace Kratos::Python
{

void  AddOperationsToPython(pybind11::module& m);

}  // namespace Kratos::Python
