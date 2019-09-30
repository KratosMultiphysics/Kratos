//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//


// System includes

// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/rans_variable_utils.h"

namespace Kratos {
namespace Python {

void AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<RansVariableUtils, VariableUtils>(m, "RansVariableUtils")
        .def(py::init<>())
        .def("ClipScalarVariable", &RansVariableUtils::ClipScalarVariable)
        .def("GetMinimumScalarValue", &RansVariableUtils::GetMinimumScalarValue)
        .def("GetMaximumScalarValue", &RansVariableUtils::GetMaximumScalarValue)
        .def("GetFlaggedMinimumScalarValue", &RansVariableUtils::GetFlaggedMinimumScalarValue)
        .def("GetFlaggedMaximumScalarValue", &RansVariableUtils::GetFlaggedMaximumScalarValue)
        .def("CopyNodalSolutionStepVariablesList",
             &RansVariableUtils::CopyNodalSolutionStepVariablesList);
}

} // namespace Python.
} // Namespace Kratos
