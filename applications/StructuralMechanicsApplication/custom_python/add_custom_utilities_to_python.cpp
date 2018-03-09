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
#include "custom_utilities/adjoint_utilities/finite_differences_utilities.h" //M.Fusseder TODO: maybe remove this (used only for controlling results)

#include "custom_processes/output_primal_solution_process.h"
#include "custom_processes/input_primal_solution_process.h"
#include "custom_processes/adjoint_processes/replace_elements_and_conditions_for_adjoint_problem_process.h"

//Response Functions
#include "custom_utilities/adjoint_utilities/structural_response_function.h"
#include "custom_utilities/adjoint_utilities/local_stress_response_function.h"
#include "custom_utilities/adjoint_utilities/nodal_displacement_response_function.h"
#include "custom_utilities/adjoint_utilities/strain_energy_response_function.h"


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
    class_<StructuralResponseFunction, boost::noncopyable>("StructuralResponseFunction", init<ModelPart&, Parameters&>())
        .def("Initialize", &StructuralResponseFunction::Initialize)
        .def("FinalizeSolutionStep", &StructuralResponseFunction::FinalizeSolutionStep)
        .def("CalculateValue", &StructuralResponseFunction::CalculateValue);

    class_<LocalStressResponseFunction, bases<StructuralResponseFunction>, boost::noncopyable>
      ("LocalStressResponseFunction", init<ModelPart&, Parameters&>());

    class_<NodalDisplacementResponseFunction, bases<StructuralResponseFunction>, boost::noncopyable>
      ("NodalDisplacementResponseFunction", init<ModelPart&, Parameters&>());

    class_<StrainEnergyResponseFunction, bases<StructuralResponseFunction>, boost::noncopyable>
      ("StrainEnergyResponseFunction", init<ModelPart&, Parameters&>());

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

