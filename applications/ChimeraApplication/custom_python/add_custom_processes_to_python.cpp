//
// ==============================================================================
//  ChimeraApplication
//
//  License:         BSD License
//                   license: ChimeraApplication/license.txt
//
//  Main authors:    Aditya Ghantasala, https://github.com/adityaghantasala
//                   Navaneeth K Narayanan
//
// ==============================================================================

// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "custom_python/add_custom_processes_to_python.h"

//#include "custom_processes/spalart_allmaras_turbulence_model_for_chimera.h"
#include "custom_processes/apply_chimera_process_monolithic.h"
#include "custom_processes/apply_chimera_process_fractional_step.h"
#include "custom_processes/rotate_region_process.h"
namespace Kratos
{

namespace Python
{

void AddCustomProcessesToPython(pybind11::module &m)
{

    namespace py = pybind11;

    py::class_<ApplyChimeraProcessMonolithic<2>, ApplyChimeraProcessMonolithic<2>::Pointer, Process>(m, "ApplyChimeraProcessMonolithic2d")
        .def(py::init<ModelPart &, Parameters>());

    py::class_<ApplyChimeraProcessMonolithic<3>, ApplyChimeraProcessMonolithic<3>::Pointer, Process>(m, "ApplyChimeraProcessMonolithic3d")
        .def(py::init<ModelPart &, Parameters>());

    py::class_<ApplyChimeraProcessFractionalStep<2>, ApplyChimeraProcessFractionalStep<2>::Pointer, Process>(m, "ApplyChimeraProcessFractionalStep2d")
        .def(py::init<ModelPart &, Parameters>());

    py::class_<ApplyChimeraProcessFractionalStep<3>, ApplyChimeraProcessFractionalStep<3>::Pointer, Process>(m, "ApplyChimeraProcessFractionalStep3d")
        .def(py::init<ModelPart &, Parameters>());

    py::class_<RotateRegionProcess, RotateRegionProcess::Pointer, Process>(m, "RotateRegionProcess")
        .def(py::init<ModelPart &, Parameters>())
        .def("SetCentreOfRotation", &RotateRegionProcess::SetCentreOfRotation)
        .def("ChangeAngularVelocity", &RotateRegionProcess::ChangeAngularVelocity);
}

} // namespace Python.

} // Namespace Kratos
