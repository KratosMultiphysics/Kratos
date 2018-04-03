// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_response_functions_to_python.h"

//Utilities
#include "custom_response_functions/custom_utilities/finite_differences_utilities.h" //M.Fusseder TODO: maybe remove this (used only for controlling results)

#include "custom_response_functions/adjoint_processes/replace_elements_and_conditions_for_adjoint_problem_process.h"

//Response Functions
#include "custom_response_functions/response_utilities/adjoint_structural_response_function.h"
#include "custom_response_functions/response_utilities/adjoint_local_stress_response_function.h"
#include "custom_response_functions/response_utilities/adjoint_nodal_displacement_response_function.h"
#include "custom_response_functions/response_utilities/adjoint_strain_energy_response_function.h"

#include "custom_response_functions/response_utilities/response_function.h"
#include "custom_response_functions/response_utilities/strain_energy_response_function.h"
#include "custom_response_functions/response_utilities/mass_response_function.h"
#include "custom_response_functions/response_utilities/eigenfrequency_response_function.h"
#include "custom_response_functions/response_utilities/eigenfrequency_response_function_lin_scal.h"


namespace Kratos
{
namespace Python
{

void  AddCustomResponseFunctionsToPython()
{
    using namespace boost::python;

    /// Processes
    class_<ReplaceElementsAndConditionsForAdjointProblemProcess , bases<Process>, boost::noncopyable >("ReplaceElementsAndConditionsForAdjointProblemProcess",
            init<ModelPart&, Parameters>());

    //Response Functions
    class_<AdjointStructuralResponseFunction, boost::noncopyable>
      ("AdjointStructuralResponseFunction", init<ModelPart&, Parameters&>())
      .def("Initialize", &AdjointStructuralResponseFunction::Initialize)
      .def("FinalizeSolutionStep", &AdjointStructuralResponseFunction::FinalizeSolutionStep)
      .def("CalculateValue", &AdjointStructuralResponseFunction::CalculateValue);

    class_<AdjointLocalStressResponseFunction, bases<AdjointStructuralResponseFunction>, boost::noncopyable>
      ("AdjointLocalStressResponseFunction", init<ModelPart&, Parameters&>());

    class_<AdjointNodalDisplacementResponseFunction, bases<AdjointStructuralResponseFunction>, boost::noncopyable>
      ("AdjointNodalDisplacementResponseFunction", init<ModelPart&, Parameters&>());

    class_<AdjointStrainEnergyResponseFunction, bases<AdjointStructuralResponseFunction>, boost::noncopyable>
      ("AdjointStrainEnergyResponseFunction", init<ModelPart&, Parameters&>());

    class_<ResponseFunction, boost::noncopyable >
      ("ResponseFunction", no_init)
      .def("Initialize", &ResponseFunction::Initialize)
      .def("CalculateValue", &ResponseFunction::CalculateValue)
      .def("CalculateGradient", &ResponseFunction::CalculateGradient);

    class_<StrainEnergyResponseFunction, bases<ResponseFunction>, boost::noncopyable >
      ("StrainEnergyResponseFunction", init<ModelPart&, Parameters>());

    class_<MassResponseFunction, bases<ResponseFunction>, boost::noncopyable >
      ("MassResponseFunction", init<ModelPart&, Parameters>());

    class_<EigenfrequencyResponseFunction, bases<ResponseFunction>, boost::noncopyable >
      ("EigenfrequencyResponseFunction", init<ModelPart&, Parameters&>());

    class_<EigenfrequencyResponseFunctionLinScal, bases<ResponseFunction>, boost::noncopyable >
      ("EigenfrequencyResponseFunctionLinScal", init<ModelPart&, Parameters&>());

    //For global finite differences
    class_<FiniteDifferencesUtilities, boost::noncopyable>("FiniteDifferencesUtilities", init< >())
        .def("SetDesignVariable", &FiniteDifferencesUtilities::SetDesignVariable)
        .def("GetDesignVariable", &FiniteDifferencesUtilities::GetDesignVariable)
        .def("SetDerivedObject", &FiniteDifferencesUtilities::SetDerivedObject)
        .def("GetDerivedObject", &FiniteDifferencesUtilities::GetDerivedObject)
        .def("DisturbElementDesignVariable", &FiniteDifferencesUtilities::DisturbElementDesignVariable)
        .def("UndisturbElementDesignVariable", &FiniteDifferencesUtilities::UndisturbElementDesignVariable)
        .def("GetStressResultantBeam", &FiniteDifferencesUtilities::GetStressResultantBeam)
        .def("GetStressResultantShell", &FiniteDifferencesUtilities::GetStressResultantShell)
        .def("GetNodalDisplacement", &FiniteDifferencesUtilities::GetNodalDisplacement)
        .def("GetStrainEnergy", &FiniteDifferencesUtilities::GetStrainEnergy);
}

}  // namespace Python.

} // Namespace Kratos

