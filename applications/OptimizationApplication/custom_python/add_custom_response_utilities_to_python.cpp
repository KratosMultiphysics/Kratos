//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

// System includes

// External includes
#include <pybind11/stl.h>

// Project includes

// Application includes
#include "custom_utilities/response/mass_response_utils.h"
#include "custom_utilities/response/linear_strain_energy_response_utils.h"

// Include base h
#include "add_custom_response_utilities_to_python.h"

namespace Kratos {
namespace Python {

void  AddCustomResponseUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    m.def_submodule("MassResponseUtils")
        .def("Check", &MassResponseUtils::Check)
        .def("CalculateValue", &MassResponseUtils::CalculateValue)
        .def("CalculateGradient", &MassResponseUtils::CalculateGradient, py::arg("evaluated_model_part"), py::arg("gradient_variable_model_parts_map"))
        ;

    m.def_submodule("LinearStrainEnergyResponseUtils")
        .def("CalculateValue", &LinearStrainEnergyResponseUtils::CalculateValue)
        .def("CalculateGradient", &LinearStrainEnergyResponseUtils::CalculateGradient, py::arg("analysis_model_part"), py::arg("evaluated_model_part"), py::arg("gradient_variable_model_parts_map"), py::arg("perturbation_size"))
        ;

}

}  // namespace Python.
} // Namespace Kratos

