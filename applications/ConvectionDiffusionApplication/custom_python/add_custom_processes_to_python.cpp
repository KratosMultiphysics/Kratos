// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Aniol Sala Pascual
//

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "processes/process.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "includes/kratos_parameters.h"

// Response Functions
#include "custom_processes/manufactured_body_force_process.h"
#include "custom_processes/microfluidic_tube_process.h"


namespace Kratos {
namespace Python {

void  AddCustomProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    // Processes
    py::class_<ManufacturedBodyForceProcess, ManufacturedBodyForceProcess::Pointer, Process>
    (m, "ManufacturedBodyForceProcess")
    .def(py::init< ModelPart&>())
    .def(py::init< ModelPart&, Parameters& >())
    ;

    py::class_<MicrofluidicTubeProcess, MicrofluidicTubeProcess::Pointer, Process>
    (m, "MicrofluidicTubeProcess")
    .def(py::init< ModelPart&>())
    .def(py::init< ModelPart&, Parameters& >())
    ;
}

}  // namespace Python.
} // Namespace Kratos