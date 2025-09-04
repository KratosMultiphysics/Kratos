
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
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"

// Legacy Extensions
#include "legacy/python/add_processes_to_python.h"

namespace Kratos::Legacy::Python
{

namespace py = pybind11;

void AddLegacyToPython(py::module& m)
{
    AddProcessesToPython(m);
}

}  // namespace Kratos::Legacy::Python.

