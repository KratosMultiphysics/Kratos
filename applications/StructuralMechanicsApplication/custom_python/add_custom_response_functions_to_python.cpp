// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:        BSD License
//	                license: structural_mechanics_application/license.txt
//
//  Main authors:    Armin Geiser
//

// System includes

// External includes
#include <pybind11/stl.h>

// Project includes
#include "custom_python/add_custom_response_functions_to_python.h"

// Processes
#include "custom_response_functions/adjoint_processes/replace_multiple_elements_and_conditions_process.h"

// Response Functions
#include "custom_response_functions/response_utilities/strain_energy_response_function_utility.h"
#include "custom_response_functions/response_utilities/mass_response_function_utility.h"
#include "custom_response_functions/response_utilities/eigenfrequency_response_function_utility.h"

#include "response_functions/adjoint_response_function.h"
#include "custom_response_functions/response_utilities/adjoint_local_stress_response_function.h"
#include "custom_response_functions/response_utilities/adjoint_max_stress_response_function.h"
#include "custom_response_functions/response_utilities/adjoint_nodal_displacement_response_function.h"
#include "custom_response_functions/response_utilities/adjoint_linear_strain_energy_response_function.h"
#include "custom_response_functions/response_utilities/adjoint_nodal_reaction_response_function.h"

namespace Kratos {
namespace Python {

void  AddCustomResponseFunctionUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    // Response Functions
    py::class_<StrainEnergyResponseFunctionUtility, StrainEnergyResponseFunctionUtility::Pointer >
        (m, "StrainEnergyResponseFunctionUtility")
        .def(py::init<ModelPart&, Parameters>())
        .def("Initialize", &StrainEnergyResponseFunctionUtility::Initialize)
        .def("CalculateValue", &StrainEnergyResponseFunctionUtility::CalculateValue)
        .def("CalculateGradient", &StrainEnergyResponseFunctionUtility::CalculateGradient);

    py::class_<MassResponseFunctionUtility, MassResponseFunctionUtility::Pointer >
        (m, "MassResponseFunctionUtility")
        .def(py::init<ModelPart&, Parameters>())
        .def("Initialize", &MassResponseFunctionUtility::Initialize)
        .def("CalculateValue", &MassResponseFunctionUtility::CalculateValue)
        .def("CalculateGradient", &MassResponseFunctionUtility::CalculateGradient);

    py::class_<EigenfrequencyResponseFunctionUtility, EigenfrequencyResponseFunctionUtility::Pointer >
        (m, "EigenfrequencyResponseFunctionUtility")
        .def(py::init<ModelPart&, Parameters>())
        .def("Initialize", &EigenfrequencyResponseFunctionUtility::Initialize)
        .def("CalculateValue", &EigenfrequencyResponseFunctionUtility::CalculateValue)
        .def("CalculateGradient", &EigenfrequencyResponseFunctionUtility::CalculateGradient);

    // Processes
    py::class_<ReplaceMultipleElementsAndConditionsProcess, ReplaceMultipleElementsAndConditionsProcess::Pointer , Process>
        (m, "ReplaceMultipleElementsAndConditionsProcess")
        .def(py::init<ModelPart&, Parameters>());

    // Response Functions
    py::class_<AdjointLocalStressResponseFunction, AdjointLocalStressResponseFunction::Pointer, AdjointResponseFunction>
        (m, "AdjointLocalStressResponseFunction")
        .def(py::init<ModelPart&, Parameters>());

    py::class_<AdjointMaxStressResponseFunction, AdjointMaxStressResponseFunction::Pointer, AdjointResponseFunction>
        (m, "AdjointMaxStressResponseFunction")
        .def(py::init<ModelPart&, Parameters>());

    py::class_<AdjointNodalDisplacementResponseFunction, AdjointNodalDisplacementResponseFunction::Pointer, AdjointResponseFunction>
        (m, "AdjointNodalDisplacementResponseFunction")
        .def(py::init<ModelPart&, Parameters>());

    py::class_<AdjointLinearStrainEnergyResponseFunction, AdjointLinearStrainEnergyResponseFunction::Pointer, AdjointResponseFunction>
        (m, "AdjointLinearStrainEnergyResponseFunction")
        .def(py::init<ModelPart&, Parameters>());

    py::class_<AdjointNodalReactionResponseFunction, AdjointNodalReactionResponseFunction::Pointer, AdjointResponseFunction>
        (m, "AdjointNodalReactionResponseFunction")
        .def(py::init<ModelPart&, Parameters>());
}

}  // namespace Python.
} // Namespace Kratos

