//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Tosi
//


// System includes


// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define_python.h"
#include "custom_python/add_custom_strategies_to_python.h"

//statistics
#include "custom_statistics/power_sums_statistics.h"


namespace Kratos {
namespace Python {

void  AddCustomStatisticsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<PowerSumsStatistics, PowerSumsStatistics::Pointer, Process >
        (m, "PowerSumsStatistics")
        .def(py::init<ModelPart&>())
        .def(py::init<ModelPart&, Parameters>())
    ;

}

} // namespace Python.
} // Namespace Kratos
