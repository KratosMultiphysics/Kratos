//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         SystemIdentificationApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//                   Ihar Antonau
//                   Fabian Meister
//

// System includes

// External includes
#include <pybind11/stl.h>

// Project includes
#include "includes/define.h"

// Application includes
#include "custom_utilities/control_utils.h"
#include "custom_utilities/smooth_clamper.h"

// Include base h
#include "custom_python/add_custom_utilities_to_python.h"

namespace Kratos::Python {

template<class TContainerType>
void AddSmoothClamper(
    pybind11::module& m,
    const std::string& rName)
{
    namespace py = pybind11;

    using smooth_clamper_type = SmoothClamper<TContainerType>;
    py::class_<smooth_clamper_type, typename smooth_clamper_type::Pointer>(m, (rName + "SmoothClamper").c_str())
        .def(py::init<const double, const double>(), py::arg("min"), py::arg("max"))
        .def("ProjectForward", &smooth_clamper_type::ProjectForward, py::arg("x_expression"))
        .def("CalculateForwardProjectionGradient", &smooth_clamper_type::CalculateForwardProjectionGradient, py::arg("x_expression"))
        .def("ProjectBackward", &smooth_clamper_type::ProjectBackward, py::arg("y_expression"))
        ;
}

void AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto control_utils = m.def_submodule("ControlUtils");
    control_utils.def("AssignEquivalentProperties", &ControlUtils::AssignEquivalentProperties<ModelPart::ConditionsContainerType>, py::arg("source_conditions"), py::arg("destination_conditions"));
    control_utils.def("AssignEquivalentProperties", &ControlUtils::AssignEquivalentProperties<ModelPart::ElementsContainerType>, py::arg("source_elements"), py::arg("destination_elements"));
    control_utils.def("ClipContainerExpression", &ControlUtils::ClipContainerExpression<ModelPart::NodesContainerType>, py::arg("nodal_expression"), py::arg("min"), py::arg("max"));
    control_utils.def("ClipContainerExpression", &ControlUtils::ClipContainerExpression<ModelPart::ConditionsContainerType>, py::arg("condition_expression"), py::arg("min"), py::arg("max"));
    control_utils.def("ClipContainerExpression", &ControlUtils::ClipContainerExpression<ModelPart::ElementsContainerType>, py::arg("element_expression"), py::arg("min"), py::arg("max"));

    AddSmoothClamper<ModelPart::NodesContainerType>(m, "Node");
    AddSmoothClamper<ModelPart::ConditionsContainerType>(m, "Condition");
    AddSmoothClamper<ModelPart::ElementsContainerType>(m, "Element");
}

} // namespace Kratos::Python
