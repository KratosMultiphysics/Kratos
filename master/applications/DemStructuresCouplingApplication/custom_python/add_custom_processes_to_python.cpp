//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//


// External includes

// Project includes
#include "custom_python/add_custom_processes_to_python.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "includes/kratos_parameters.h"

#include "custom_processes/control_module_process.hpp"


namespace Kratos
{

namespace Python
{

void  AddCustomProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<ControlModuleProcess, ControlModuleProcess::Pointer, Process>
    (m, "ControlModuleProcess")
    .def( py::init< ModelPart&, Parameters>());
}

}  // namespace Python.
} // Namespace Kratos
