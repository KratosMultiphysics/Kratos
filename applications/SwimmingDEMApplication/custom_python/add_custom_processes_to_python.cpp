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
//                   Guillermo Casas
//                   Joaquin Gonzalez-Usua
//

// External includes

// Project includes
#include "includes/model_part.h"
#include "processes/process.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "includes/kratos_parameters.h"

#include "custom_processes/apply_rigid_rotation_process.hpp"
#include "custom_processes/bump_transient_porosity_solution_body_force_process.h"
#include "custom_processes/porosity_solution_and_sinusoidal_body_force_process.h"
#include "custom_processes/porosity_solution_and_body_force_process.h"
#include "custom_processes/sinusoidal_porosity_solution_transient_body_force_process.h"
#include "custom_processes/sinusoidal_porosity_solution_and_body_force_process.h"
#include "custom_processes/hyperbolic_tangential_porosity_solution_and_body_force_process.h"
#include "custom_processes/hyperbolic_tangential_porosity_solution_transient_body_force_process.h"
#include "custom_processes/porosity_solution_transient_body_force_process.h"

namespace Kratos
{

namespace Python
{

void  AddCustomProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<ApplyRigidRotationProcess, ApplyRigidRotationProcess::Pointer, Process>
    (m, "ApplyRigidRotationProcess")
    .def( py::init< ModelPart&, Parameters&>());

    py::class_<BumpTransientPorositySolutionBodyForceProcess, BumpTransientPorositySolutionBodyForceProcess::Pointer, Process>
    (m, "BumpTransientPorositySolutionBodyForceProcess")
    .def(py::init< ModelPart&>())
    .def(py::init< ModelPart&, Parameters& >())
    ;

    py::class_<PorositySolutionAndBodyForceProcess, PorositySolutionAndBodyForceProcess::Pointer, Process>
    (m, "PorositySolutionAndBodyForceProcess")
    .def(py::init< ModelPart&>())
    .def(py::init< ModelPart&, Parameters& >())
    ;

        py::class_<PorositySolutionAndSinusoidalBodyForceProcess, PorositySolutionAndSinusoidalBodyForceProcess::Pointer, Process>
    (m, "PorositySolutionAndSinusoidalBodyForceProcess")
    .def(py::init< ModelPart&>())
    .def(py::init< ModelPart&, Parameters& >())
    ;

    py::class_<PorositySolutionTransientBodyForceProcess, PorositySolutionTransientBodyForceProcess::Pointer, Process>
    (m, "PorositySolutionTransientBodyForceProcess")
    .def(py::init< ModelPart&>())
    .def(py::init< ModelPart&, Parameters& >())
    ;

    py::class_<SinusoidalPorositySolutionTransientBodyForceProcess, SinusoidalPorositySolutionTransientBodyForceProcess::Pointer, Process>
    (m, "SinusoidalPorositySolutionTransientBodyForceProcess")
    .def(py::init< ModelPart&>())
    .def(py::init< ModelPart&, Parameters& >())
    ;

    py::class_<SinusoidalPorositySolutionAndBodyForceProcess, SinusoidalPorositySolutionAndBodyForceProcess::Pointer, Process>
    (m, "SinusoidalPorositySolutionAndBodyForceProcess")
    .def(py::init< ModelPart&>())
    .def(py::init< ModelPart&, Parameters& >())
    ;

    py::class_<HyperbolicTangentialPorositySolutionAndBodyForceProcess, HyperbolicTangentialPorositySolutionAndBodyForceProcess::Pointer, Process>
    (m, "HyperbolicTangentialPorositySolutionAndBodyForceProcess")
    .def(py::init< ModelPart&>())
    .def(py::init< ModelPart&, Parameters& >())
    ;

    py::class_<HyperbolicTangentialPorositySolutionTransientBodyForceProcess, HyperbolicTangentialPorositySolutionTransientBodyForceProcess::Pointer, Process>
    (m, "HyperbolicTangentialPorositySolutionTransientBodyForceProcess")
    .def(py::init< ModelPart&>())
    .def(py::init< ModelPart&, Parameters& >())
    ;
}

}  // namespace Python.
} // Namespace Kratos
