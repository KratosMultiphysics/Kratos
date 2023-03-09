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
#include "custom_utilities/response/response_utils.h"
#include "custom_utilities/response/mass_response_utils.h"
#include "custom_utilities/response/linear_strain_energy_response_utils.h"

// Include base h
#include "add_custom_response_utilities_to_python.h"

namespace Kratos {
namespace Python {

void  AddCustomResponseUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    m.def("GetSensitivityModelPartForAdjointSensitivities", &ResponseUtils::GetSensitivityModelPartForAdjointSensitivities,
        py::arg("sensitivity_model_parts_list"),
        py::arg("analysis_model_part"),
        py::arg("are_sensitivity_entity_parents_considered"),
        py::arg("are_sensitivity_entites_considered"),
        py::arg("force_find_sensitivity_entities_in_analysis_model_part") = false);
    m.def("GetSensitivityModelPartForDirectSensitivities", &ResponseUtils::GetSensitivityModelPartForDirectSensitivities,
        py::arg("sensitivity_model_parts_list"),
        py::arg("evaluated_model_parts_list"),
        py::arg("are_nodes_considered"),
        py::arg("are_conditions_considered"),
        py::arg("are_elements_considered"));

    py::class_<MassResponseUtils >(m, "MassResponseUtils")
        .def_static("Check", &MassResponseUtils::Check)
        .def_static("CalculateValue", &MassResponseUtils::CalculateValue)
        .def_static("CalculateSensitivity", &MassResponseUtils::CalculateSensitivity)
        ;

    py::class_<LinearStrainEnergyResponseUtils >(m, "LinearStrainEnergyResponseUtils")
        .def_static("CalculateValue", &LinearStrainEnergyResponseUtils::CalculateValue)
        .def_static("CalculateSensitivity", &LinearStrainEnergyResponseUtils::CalculateSensitivity)
        ;

}

}  // namespace Python.
} // Namespace Kratos

