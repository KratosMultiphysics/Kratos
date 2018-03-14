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
#include "custom_python/add_custom_utilities_to_python.h"

//Utilities
#include "custom_utilities/sprism_neighbours.hpp"
#include "custom_response_functions/custom_utilities/finite_differences_utilities.h" //M.Fusseder TODO: maybe remove this (used only for controlling results)

#include "custom_response_functions/adjoint_processes/output_primal_solution_process.h"
#include "custom_response_functions/adjoint_processes/input_primal_solution_process.h"
#include "custom_response_functions/adjoint_processes/replace_elements_and_conditions_for_adjoint_problem_process.h"

//Response Functions
#include "custom_response_functions/response_utilities/adjoint_structural_response_function.h"
#include "custom_response_functions/response_utilities/adjoint_local_stress_response_function.h"
#include "custom_response_functions/response_utilities/adjoint_nodal_displacement_response_function.h"
#include "custom_response_functions/response_utilities/adjoint_strain_energy_response_function.h"

#include "custom_response_functions/response_utilities/strain_energy_response_function.h"
#include "custom_response_functions/response_utilities/mass_response_function.h"
#include "custom_response_functions/response_utilities/eigenfrequency_response_function.h"
#include "custom_response_functions/response_utilities/eigenfrequency_response_function_lin_scal.h"


namespace Kratos
{
namespace Python
{

void  AddCustomUtilitiesToPython()
{
    using namespace boost::python;

//     typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
//     typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
//     typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

    class_<SprismNeighbours>("SprismNeighbours", init<ModelPart&>())
    .def("Execute",&SprismNeighbours::Execute)
    .def("ClearNeighbours",&SprismNeighbours::ClearNeighbours)
    ;

    /// Processes
    class_< OutputPrimalSolutionProcess, bases<Process> >
    ("OutputPrimalSolutionProcess", init<ModelPart&, Parameters&>());

    class_< InputPrimalSolutionProcess, bases<Process> >
    ("InputPrimalSolutionProcess", init<ModelPart&, Parameters&>());

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

    class_<StrainEnergyResponseFunction, boost::noncopyable >
      ("StrainEnergyResponseFunction", init<ModelPart&, Parameters>())
      .def("Initialize", &StrainEnergyResponseFunction::Initialize)
      .def("CalculateValue", &StrainEnergyResponseFunction::CalculateValue)
      .def("CalculateGradient", &StrainEnergyResponseFunction::CalculateGradient)
      .def("GetValue", &StrainEnergyResponseFunction::GetValue)
      .def("GetGradient", &StrainEnergyResponseFunction::GetGradient);

    class_<MassResponseFunction, boost::noncopyable >
      ("MassResponseFunction", init<ModelPart&, Parameters>())
      .def("Initialize", &MassResponseFunction::Initialize)
      .def("CalculateValue", &MassResponseFunction::CalculateValue)
      .def("CalculateGradient", &MassResponseFunction::CalculateGradient)
      .def("GetValue", &MassResponseFunction::GetValue)
      .def("GetGradient", &MassResponseFunction::GetGradient);

    class_<EigenfrequencyResponseFunction, boost::noncopyable >
      ("EigenfrequencyResponseFunction", init<ModelPart&, Parameters&>())
      .def("Initialize", &EigenfrequencyResponseFunction::Initialize)
      .def("CalculateValue", &EigenfrequencyResponseFunction::CalculateValue)
      .def("CalculateGradient", &EigenfrequencyResponseFunction::CalculateGradient)
      .def("GetValue", &EigenfrequencyResponseFunction::GetValue)
      .def("GetGradient", &EigenfrequencyResponseFunction::GetGradient);

    class_<EigenfrequencyResponseFunctionLinScal, boost::noncopyable >
      ("EigenfrequencyResponseFunctionLinScal", init<ModelPart&, Parameters&>())
      .def("Initialize", &EigenfrequencyResponseFunctionLinScal::Initialize)
      .def("CalculateValue", &EigenfrequencyResponseFunctionLinScal::CalculateValue)
      .def("CalculateGradient", &EigenfrequencyResponseFunctionLinScal::CalculateGradient)
      .def("GetValue", &EigenfrequencyResponseFunctionLinScal::GetValue)
      .def("GetGradient", &EigenfrequencyResponseFunctionLinScal::GetGradient);

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

