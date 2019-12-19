//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Angel Celigueta
//
//

// External includes

// Project includes
#include "custom_python/add_custom_processes_to_python.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "includes/kratos_parameters.h"

// Processes
#include "custom_processes/apply_kinematic_constraints_process.hpp"
#include "custom_processes/control_module_2d_process.hpp"

namespace Kratos
{

namespace Python
{

void  AddCustomProcessesToPython(pybind11::module& m)
{

    namespace py = pybind11;

    // Apply table values
    py::class_<ApplyKinematicConstraintsProcess, ApplyKinematicConstraintsProcess::Pointer, Process>
    (m, "ApplyKinematicConstraintsProcess")
    .def(py::init < ModelPart&, Parameters>());

    py::class_<ControlModule2DProcess, ControlModule2DProcess::Pointer, Process>
    (m, "ControlModule2DProcess")
    .def( py::init< ModelPart&, Parameters>());

}

}  // namespace Python.
} // Namespace Kratos

