//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//

#pragma once

// System includes

// External includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// Project includes

namespace Kratos::Python
{

void AddPropertiesToPython(pybind11::module& m);

}  // namespace Kratos::Python.
