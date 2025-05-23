// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Mate Kelemen
//

#pragma once

// External includes
#include <pybind11/pybind11.h>

namespace Kratos::Python
{
    void AddCustomConstraintsToPython(pybind11::module& rModule);
}  // namespace Kratos::Python.
