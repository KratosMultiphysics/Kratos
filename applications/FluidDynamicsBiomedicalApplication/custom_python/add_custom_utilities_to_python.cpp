//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Eduardo Soudah
//
//


// System includes

// External includes

// Project includes
#include "custom_python/add_custom_utilities_to_python.h"

// Application includes
#include "custom_utilities/wss_statistics_utilities.h"


namespace Kratos
{

namespace Python
{

void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    // WSS statistics utilities
    py::class_<WssStatisticsUtilities>(m,"WssStatisticsUtilities")
        .def(py::init<>())
        .def("CalculateWSS", &WssStatisticsUtilities::CalculateWSS)
        .def("CalculateTWSS", &WssStatisticsUtilities::CalculateTWSS)
        .def("CalculateOSI", &WssStatisticsUtilities::CalculateOSI)
        ;

}

}  // namespace Python.

} // Namespace Kratos
