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
#include "custom_processes/apply_velocity_constraints_process.hpp"
#include "custom_processes/apply_angular_velocity_constraints_process.hpp"
#include "custom_processes/apply_forces_process.hpp"
#include "custom_processes/apply_moments_process.hpp"
#include "custom_processes/control_module_2d_process.hpp"
#include "custom_processes/automatic_dt_process.hpp"
#include "custom_processes/apply_particle_injection_process.hpp"

namespace Kratos
{

namespace Python
{

void  AddCustomProcessesToPython(pybind11::module& m)
{

    namespace py = pybind11;

    // Apply table values
    py::class_<ApplyVelocityConstraintsProcess, ApplyVelocityConstraintsProcess::Pointer, Process>
    (m, "ApplyVelocityConstraintsProcess")
    .def(py::init < ModelPart&, Parameters>());

    py::class_<ApplyAngularVelocityConstraintsProcess, ApplyAngularVelocityConstraintsProcess::Pointer, Process>
    (m, "ApplyAngularVelocityConstraintsProcess")
    .def(py::init < ModelPart&, Parameters>());

    py::class_<ApplyMomentsProcess, ApplyMomentsProcess::Pointer, Process>
    (m, "ApplyMomentsProcess")
    .def(py::init < ModelPart&, Parameters >());

    py::class_<ApplyForcesProcess, ApplyForcesProcess::Pointer, Process>
    (m, "ApplyForcesProcess")
    .def(py::init < ModelPart&, Parameters >());

    py::class_<ControlModule2DProcess, ControlModule2DProcess::Pointer, Process>
    (m, "ControlModule2DProcess")
    .def( py::init< ModelPart&, Parameters>());

    py::class_<AutomaticDTProcess, AutomaticDTProcess::Pointer, Process>
    (m, "AutomaticDTProcess")
    .def( py::init< ModelPart&, Parameters>());

    py::class_<ApplyParticleInjectionProcess, ApplyParticleInjectionProcess::Pointer, Process>
    (m, "ApplyParticleInjectionProcess")
    .def(py::init < ModelPart&, Parameters >());

}

}  // namespace Python.
} // Namespace Kratos

