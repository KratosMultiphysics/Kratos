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
//

// External includes

// Project includes
#include "includes/model_part.h"
#include "processes/process.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "includes/kratos_parameters.h"

#include "custom_processes/apply_rigid_rotation_process.hpp"
#include "custom_processes/transient_spatial_dependant_porosity_solution_body_force_process.h"
#include "custom_processes/spatial_dependant_porosity_solution_body_force_process.h"

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

    py::class_<TransientSpatialDependantPorositySolutionBodyForceProcess, TransientSpatialDependantPorositySolutionBodyForceProcess::Pointer, Process>
    (m, "TransientSpatialDependantPorositySolutionBodyForceProcess")
    .def(py::init< ModelPart&, const double, const double, const double, const double, const double, const double, const double, const double>())
    .def(py::init< ModelPart&, Parameters& >())
    ;

    py::class_<SpatialDependantPorositySolutionBodyForceProcess, SpatialDependantPorositySolutionBodyForceProcess::Pointer, Process>
    (m, "SpatialDependantPorositySolutionBodyForceProcess")
    .def(py::init< ModelPart&, const double, const double, const double, const double, const double, const double>())
    .def(py::init< ModelPart&, Parameters& >())
    ;
}

}  // namespace Python.
} // Namespace Kratos
