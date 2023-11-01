//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//                   Suneth Warnakulasuriya
//

#pragma once

// System includes

// External includes
#include "pybind11/pybind11.h"

// Project includes

namespace Kratos {
namespace Python {

void  AddCustomOperationsToPython(pybind11::module& m);

}  // namespace Python.
}  // namespace Kratos.