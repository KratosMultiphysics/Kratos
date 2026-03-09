//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         SystemIdentificationApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//                   Ihar Antonau
//                   Fabian Meister
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define_python.h"

namespace Kratos::Python
{

void AddCustomUtilitiesToPython(pybind11::module& m);

} // namespace Kratos::Python
