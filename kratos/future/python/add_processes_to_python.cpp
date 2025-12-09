
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
#include "future/python/add_processes_to_python.h"
#include "future/processes/future_process.h"

namespace Kratos::Future::Python
{

namespace py = pybind11;

void AddProcessesToPython(py::module& m)
{
    // Pack everything under the "Future" submodule
    py::class_<Future::Process, Future::Process::Pointer, Flags>(m,"Process")
        .def(py::init<>())
        .def("Execute",&Future::Process::Execute)
        .def("Info",&Future::Process::Info)
        .def("__str__", PrintObject<Future::Process>)
    ;
}

}  // namespace Kratos::Future::Python.

