
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
#include "legacy/processes/legacy_process.h"

namespace Kratos::Legacy::Python
{

namespace py = pybind11;

void AddProcessesToPython(py::module& m)
{
    // Pack everything under the "Legacy" submodule
    py::class_<Legacy::Process, Legacy::Process::Pointer, Flags>(m,"Process")
        .def(py::init<>())
        .def("Execute",&Legacy::Process::Execute)
        .def("Info",&Legacy::Process::Info)
        .def("__str__", PrintObject<Legacy::Process>)
    ;
}

}  // namespace Kratos::Python.

