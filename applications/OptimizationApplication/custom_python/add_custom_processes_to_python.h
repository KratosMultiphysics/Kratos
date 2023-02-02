//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: HDF5Application/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

#pragma once

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"

namespace Kratos {
namespace Python {

void AddCustomProcessesToPython(pybind11::module& m);

} // namespace Python.
} // namespace Kratos.
