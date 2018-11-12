//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Jordi Cotela
//

// System includes

// External includes

// Project includes
#include "add_parallel_environment_to_python.h"
#include "includes/define_python.h"
#include "includes/parallel_environment.h"

namespace Kratos {

namespace Python {

void AddParallelEnvironmentToPython(pybind11::module &m)
{
    namespace py = pybind11;

    py::class_<ParallelEnvironment, ParallelEnvironment::Pointer>(m,"ParallelEnvironment")
    .def_property_readonly("MakeDefault", [](const ParallelEnvironment& self) { return ParallelEnvironment::MakeDefault; } )
    .def_property_readonly("DoNotMakeDefault", [](const ParallelEnvironment& self) { return ParallelEnvironment::DoNotMakeDefault; } )
    .def_static("GetInstance", &ParallelEnvironment::GetInstance, py::return_value_policy::reference)
    .def("RegisterDataCommunicator", &ParallelEnvironment::RegisterDataCommunicator)
    .def("GetDataCommunicator",&ParallelEnvironment::GetDataCommunicator, py::return_value_policy::reference_internal)
    .def("GetDefaultDataCommunicator",&ParallelEnvironment::GetDefaultDataCommunicator, py::return_value_policy::reference_internal)
    .def("SetDefaultDataCommunicator",&ParallelEnvironment::SetDefaultDataCommunicator)
    .def("HasDataCommunicator",&ParallelEnvironment::HasDataCommunicator)
    .def("__str__", PrintObject<ParallelEnvironment>);

}

} // namespace Python.

} // Namespace Kratos
