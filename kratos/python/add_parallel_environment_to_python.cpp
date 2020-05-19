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
#include "includes/parallel_environment.h"

namespace Kratos {

namespace Python {

void AddParallelEnvironmentToPython(pybind11::module &m)
{
    namespace py = pybind11;

    py::class_<ParallelEnvironment, std::unique_ptr<ParallelEnvironment, py::nodelete>>(m,"ParallelEnvironment")
    .def_property_readonly("MakeDefault", [](const ParallelEnvironment& self) { return ParallelEnvironment::MakeDefault; } )
    .def_property_readonly("DoNotMakeDefault", [](const ParallelEnvironment& self) { return ParallelEnvironment::DoNotMakeDefault; } )
    .def_static("UnregisterDataCommunicator", &ParallelEnvironment::UnregisterDataCommunicator)
    .def_static("GetDataCommunicator",&ParallelEnvironment::GetDataCommunicator, py::return_value_policy::reference)
    .def_static("GetDefaultDataCommunicator",&ParallelEnvironment::GetDefaultDataCommunicator, py::return_value_policy::reference)
    .def_static("SetDefaultDataCommunicator",&ParallelEnvironment::SetDefaultDataCommunicator)
    .def_static("GetDefaultRank",&ParallelEnvironment::GetDefaultRank)
    .def_static("GetDefaultSize",&ParallelEnvironment::GetDefaultSize)
    .def_static("HasDataCommunicator",&ParallelEnvironment::HasDataCommunicator)
    .def_static("GetDefaultDataCommunicatorName",&ParallelEnvironment::GetDefaultDataCommunicatorName)
    .def_static("Info", []() {
        std::stringstream ss;
        ParallelEnvironment::PrintInfo(ss);
        ss << std::endl;
        ParallelEnvironment::PrintData(ss);
        return ss.str();
    });
}

} // namespace Python.

} // Namespace Kratos
