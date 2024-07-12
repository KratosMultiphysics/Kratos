//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//

#pragma once

// Project includes
#include <pybind11/pybind11.h>

namespace Kratos::Python
{

void  AddLoggerToPython(pybind11::module& m);

}  // namespace Kratos::Python.
