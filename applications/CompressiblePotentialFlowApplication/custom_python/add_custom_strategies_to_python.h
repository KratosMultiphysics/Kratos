//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Marco Antonio Zuñiga Perez
//                   
//

#pragma once

// System includes


// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"


namespace Kratos::Python
{

void AddCustomStrategiesToPython(pybind11::module& m);


}  // namespace Kratos::Python.

