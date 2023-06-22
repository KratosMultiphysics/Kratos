//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Reza Najian Asl
//

// System includes

// External includes
#include <pybind11/stl.h>

// Project includes

// Application includes
#include "custom_utilities/control/sigmoidal_projection_utils.h"

// Include base h
#include "add_custom_control_utilities_to_python.h"

namespace Kratos {
namespace Python {

void  AddCustomControlUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    m.def_submodule("ControlUtils")
        .def("OneLevelForwardProjection", &SigmoidalProjectionUtils::OneLevelForwardProjection<ModelPart::NodesContainerType>, py::arg("input_expression"), py::arg("projected_expression"))
        .def("OneLevelForwardProjection", &SigmoidalProjectionUtils::OneLevelForwardProjection<ModelPart::ConditionsContainerType>, py::arg("input_expression"), py::arg("projected_expression"))
        .def("OneLevelForwardProjection", &SigmoidalProjectionUtils::OneLevelForwardProjection<ModelPart::ElementsContainerType>, py::arg("input_expression"), py::arg("projected_expression"))
        ;

}

}  // namespace Python.
} // Namespace Kratos

