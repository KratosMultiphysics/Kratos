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
//                   Suneth Warnakulasuriya
//

// System includes
#include <string>
#include <type_traits>

// External includes
#include <pybind11/stl.h>

// Project includes

// Application includes
#include "custom_utilities/control/sigmoidal_projection_utils.h"
#include "custom_utilities/control/smooth_clamper.h"

// Include base h
#include "add_custom_control_utilities_to_python.h"

namespace Kratos {
namespace Python {

void  AddCustomControlUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto module = m.def_submodule("SigmoidalProjectionUtils");
    module.def("ProjectForward", &SigmoidalProjectionUtils::ProjectForward, py::arg("input_tensor_adaptor"), py::arg("x_values"), py::arg("y_values"), py::arg("beta"), py::arg("penalty_factor"));
    module.def("ProjectBackward", &SigmoidalProjectionUtils::ProjectBackward, py::arg("input_tensor_adaptor"), py::arg("x_values"), py::arg("y_values"), py::arg("beta"), py::arg("penalty_factor"));
    module.def("CalculateForwardProjectionGradient", &SigmoidalProjectionUtils::CalculateForwardProjectionGradient, py::arg("input_tensor_adaptor"), py::arg("x_values"), py::arg("y_values"), py::arg("beta"), py::arg("penalty_factor"));

    py::class_<SmoothClamper, typename SmoothClamper::Pointer>(m, "SmoothClamper")
        .def(py::init<const double, const double>(), py::arg("min"), py::arg("max"))
        .def("ProjectForward", py::overload_cast<const TensorAdaptor<double>&>(&SmoothClamper::ProjectForward, py::const_), py::arg("x_tensor_adaptor"))
        .def("CalculateForwardProjectionGradient", py::overload_cast<const TensorAdaptor<double>&>(&SmoothClamper::CalculateForwardProjectionGradient, py::const_), py::arg("x_tensor_adaptor"))
        .def("ProjectBackward", py::overload_cast<const TensorAdaptor<double>&>(&SmoothClamper::ProjectBackward, py::const_), py::arg("y_tensor_adaptor"))
        ;
}

}  // namespace Python.
} // Namespace Kratos

