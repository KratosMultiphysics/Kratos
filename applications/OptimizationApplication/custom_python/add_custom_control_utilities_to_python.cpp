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

    m.def_submodule("SigmoidalProjectionUtils")
        .def("ProjectForward", &SigmoidalProjectionUtils::ProjectForward<ModelPart::NodesContainerType>)
        .def("ProjectForward", &SigmoidalProjectionUtils::ProjectForward<ModelPart::ConditionsContainerType>)
        .def("ProjectForward", &SigmoidalProjectionUtils::ProjectForward<ModelPart::ElementsContainerType>)
        .def("ProjectBackward", &SigmoidalProjectionUtils::ProjectBackward<ModelPart::NodesContainerType>)
        .def("ProjectBackward", &SigmoidalProjectionUtils::ProjectBackward<ModelPart::ConditionsContainerType>)
        .def("ProjectBackward", &SigmoidalProjectionUtils::ProjectBackward<ModelPart::ElementsContainerType>)
        ;

}

}  // namespace Python.
} // Namespace Kratos

