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
#include <string>
#include <type_traits>

// External includes
#include <pybind11/stl.h>

// Project includes

// Application includes
#include "custom_utilities/control/sigmoidal_projection_utils.h"

// Include base h
#include "add_custom_control_utilities_to_python.h"

namespace Kratos {
namespace Python {

namespace detail
{
template <class TContainerType>
void AddSigmoidalProjectionUtils(pybind11::module& m)
{
    namespace py = pybind11;

    std::string container_type;
    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        container_type = "nodal_expression";
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
        container_type = "condition_expression";
    } else {
        container_type = "element_expression";
    }

    m.def("ProjectForward", &SigmoidalProjectionUtils::ProjectForward<TContainerType>, py::arg(container_type.c_str()), py::arg("x_values"), py::arg("y_values"), py::arg("beta"), py::arg("penalty_factor"));
    m.def("ProjectBackward", &SigmoidalProjectionUtils::ProjectBackward<TContainerType>, py::arg(container_type.c_str()), py::arg("x_values"), py::arg("y_values"), py::arg("beta"), py::arg("penalty_factor"));
    m.def("CalculateForwardProjectionGradient", &SigmoidalProjectionUtils::CalculateForwardProjectionGradient<TContainerType>, py::arg(container_type.c_str()), py::arg("x_values"), py::arg("y_values"), py::arg("beta"), py::arg("penalty_factor"));
}
} // namespace detail

void  AddCustomControlUtilitiesToPython(pybind11::module& m)
{
    auto module = m.def_submodule("SigmoidalProjectionUtils");
    detail::AddSigmoidalProjectionUtils<ModelPart::NodesContainerType>(module);
    detail::AddSigmoidalProjectionUtils<ModelPart::ConditionsContainerType>(module);
    detail::AddSigmoidalProjectionUtils<ModelPart::ElementsContainerType>(module);
}

}  // namespace Python.
} // Namespace Kratos

