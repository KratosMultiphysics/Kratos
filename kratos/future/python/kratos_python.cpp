
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

// Future Extensions
#include "future/python/add_linear_solvers_to_python.h"
#include "future/python/add_processes_to_python.h"

namespace Kratos::Future::Python
{

namespace py = pybind11;

void AddFutureToPython(py::module& m)
{
    AddLinearSolversToPython(m);

    AddProcessesToPython(m);
}

}  // namespace Kratos::Future::Python.

