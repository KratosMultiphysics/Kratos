// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
// | | | (_)        |  \/  |         | |
// | | | |_ ___  ___| .  . | ___   __| |
// | | | | / __|/ __| |\/| |/ _ \ / _` |
// \ \_/ / \__ \ (__| |  | | (_) | (_| |
//  \___/|_|___/\___\_|  |_/\___/ \__,_|  APPLICATION
//                                      
//
//  License:         BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aniol Sala
//

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"
#include "processes/process.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_processes/vm_boussinesq_concentration_field_process.h"

namespace Kratos
{

namespace Python
{

void  AddCustomProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<VmBoussinesqConcentrationFieldProcess, VmBoussinesqConcentrationFieldProcess::Pointer, Process>
    (m, "BoussinesqConcentrationFieldProcess")
    .def(py::init< ModelPart& >())
    .def(py::init<ModelPart&, Parameters& >())
    ;
}

}  // namespace Python.
} // Namespace Kratos
