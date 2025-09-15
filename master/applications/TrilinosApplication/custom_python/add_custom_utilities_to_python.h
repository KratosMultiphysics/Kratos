//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#pragma once

// System includes
#include <pybind11/pybind11.h>

// External includes

// Project includes
#include "includes/define_python.h"

namespace Kratos::Python
{
void  AddCustomUtilitiesToPython(pybind11::module& m);
}  // namespace Kratos::Python.
