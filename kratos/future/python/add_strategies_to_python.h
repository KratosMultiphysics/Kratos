
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

#pragma once

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes

namespace Kratos::Future::Python
{
    void AddStrategiesToPython(pybind11::module &m);
}  // namespace Kratos::Python.
