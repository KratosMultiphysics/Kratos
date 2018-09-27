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
#include "custom_response_functions/adjoint_processes/replace_elements_and_conditions_for_adjoint_problem_process.h"

// Response Functions
#include "custom_response_functions/response_utilities/strain_energy_response_function_utility.h"
#include "custom_response_functions/response_utilities/mass_response_function_utility.h"
#include "custom_response_functions/response_utilities/eigenfrequency_response_function_utility.h"

#include "custom_response_functions/response_utilities/adjoint_structural_response_function.h"
#include "custom_response_functions/response_utilities/adjoint_local_stress_response_function.h"
#include "custom_response_functions/response_utilities/adjoint_nodal_displacement_response_function.h"
#include "custom_response_functions/response_utilities/adjoint_linear_strain_energy_response_function.h"

// Adjoint postprocessing
#include "custom_response_functions/response_utilities/adjoint_postprocess.h"

namespace Kratos
{
namespace Python
{

using namespace pybind11;

void  AddCustomResponseFunctionUtilitiesToPython(pybind11::module& m)
{

    // Response Functions
    class_<StrainEnergyResponseFunctionUtility, StrainEnergyResponseFunctionUtility::Pointer >
      (m, "StrainEnergyResponseFunctionUtility")
      .def(init<ModelPart&, Parameters>())
      .def("Initialize", &StrainEnergyResponseFunctionUtility::Initialize)
      .def("CalculateValue", &StrainEnergyResponseFunctionUtility::CalculateValue)
      .def("CalculateGradient", &StrainEnergyResponseFunctionUtility::CalculateGradient);

    class_<MassResponseFunctionUtility, MassResponseFunctionUtility::Pointer >
      (m, "MassResponseFunctionUtility")
      .def(init<ModelPart&, Parameters>())
      .def("Initialize", &MassResponseFunctionUtility::Initialize)
      .def("CalculateValue", &MassResponseFunctionUtility::CalculateValue)
      .def("CalculateGradient", &MassResponseFunctionUtility::CalculateGradient);

    class_<EigenfrequencyResponseFunctionUtility, EigenfrequencyResponseFunctionUtility::Pointer >
      (m, "EigenfrequencyResponseFunctionUtility")
      .def(init<ModelPart&, Parameters>())
      .def("Initialize", &EigenfrequencyResponseFunctionUtility::Initialize)
      .def("CalculateValue", &EigenfrequencyResponseFunctionUtility::CalculateValue)
      .def("CalculateGradient", &EigenfrequencyResponseFunctionUtility::CalculateGradient);

    // Processes
    class_<ReplaceElementsAndConditionsForAdjointProblemProcess, ReplaceElementsAndConditionsForAdjointProblemProcess::Pointer , Process>
      (m, "ReplaceElementsAndConditionsForAdjointProblemProcess")
      .def(init<ModelPart&>());

    // Response Functions
    class_<AdjointStructuralResponseFunction, AdjointStructuralResponseFunction::Pointer>
      (m, "AdjointStructuralResponseFunction")
      .def(init<ModelPart&, Parameters>())
      .def("Initialize", &AdjointStructuralResponseFunction::Initialize)
      .def("InitializeSolutionStep", &AdjointStructuralResponseFunction::InitializeSolutionStep)
      .def("FinalizeSolutionStep", &AdjointStructuralResponseFunction::FinalizeSolutionStep)
      .def("CalculateValue", &AdjointStructuralResponseFunction::CalculateValue)
      .def("UpdateSensitivities", &AdjointStructuralResponseFunction::UpdateSensitivities);

    class_<AdjointLocalStressResponseFunction, AdjointLocalStressResponseFunction::Pointer, AdjointStructuralResponseFunction>
      (m, "AdjointLocalStressResponseFunction")
      .def(init<ModelPart&, Parameters>());

    class_<AdjointNodalDisplacementResponseFunction, AdjointNodalDisplacementResponseFunction::Pointer, AdjointStructuralResponseFunction>
      (m, "AdjointNodalDisplacementResponseFunction")
      .def(init<ModelPart&, Parameters>());

    class_<AdjointLinearStrainEnergyResponseFunction, AdjointLinearStrainEnergyResponseFunction::Pointer, AdjointStructuralResponseFunction>
      (m, "AdjointLinearStrainEnergyResponseFunction")
      .def(init<ModelPart&, Parameters>());

    // Adjoint postprocess
    class_<AdjointPostprocess, AdjointPostprocess::Pointer>
      (m, "AdjointPostprocess")
      .def(init<ModelPart&, AdjointStructuralResponseFunction&, Parameters>())
      .def("Initialize", &AdjointPostprocess::Initialize)
      .def("InitializeSolutionStep", &AdjointPostprocess::InitializeSolutionStep)
      .def("FinalizeSolutionStep", &AdjointPostprocess::FinalizeSolutionStep)
      .def("UpdateSensitivities", &AdjointPostprocess::UpdateSensitivities);
}

}  // namespace Python.

} // Namespace Kratos

