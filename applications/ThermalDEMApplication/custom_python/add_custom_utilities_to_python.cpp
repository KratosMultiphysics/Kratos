//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics ThermalDEM Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Rafael Rangel (rrangel@cimne.upc.edu)
//

// System includes

// External includes

// Project includes
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/graph_utilities.h"
#include "custom_utilities/set_thermal_data_utilities.h"
#include "custom_utilities/tesselation_utilities.h"

namespace Kratos
{
namespace Python
{

namespace py = pybind11;

void AddCustomUtilitiesToPython(pybind11::module& m) {

  py::class_<GraphUtilities, GraphUtilities::Pointer>(m, "GraphUtilities")
    .def(py::init<>())
    .def("ExecuteInitialize", &GraphUtilities::ExecuteInitialize)
    .def("ExecuteFinalizeSolutionStep", &GraphUtilities::ExecuteFinalizeSolutionStep)
    .def("ExecuteFinalize", &GraphUtilities::ExecuteFinalize)
    ;

  py::class_<SetThermalDataUtilities, SetThermalDataUtilities::Pointer>(m, "SetThermalDataUtilities")
    .def(py::init<>())
    .def("ExecuteInitialize", &SetThermalDataUtilities::ExecuteInitialize)
    ;

  py::class_<TesselationUtilities, TesselationUtilities::Pointer>(m, "TesselationUtilities")
    .def(py::init<>())
    .def("ExecuteInitialize", &TesselationUtilities::ExecuteInitialize)
    .def("ExecuteInitializeSolutionStep", &TesselationUtilities::ExecuteInitializeSolutionStep)
    ;
}

} // namespace Python
} // namespace Kratos
