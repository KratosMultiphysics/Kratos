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
#include "includes/model_part.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

//Utilities
#include "custom_utilities/sprism_neighbours.hpp"
#include "custom_utilities/eigenvector_to_solution_step_variable_transfer_utility.hpp"
#include "custom_utilities/adjoint_utilities/finite_differences_utilities.h" //MFusseder TODO: maybe remove this only for controlling results

//Processes
#include "custom_processes/apply_multi_point_constraints_process.h"
#include "custom_processes/output_primal_solution_process.h"
#include "custom_processes/input_primal_solution_process.h"
#include "custom_processes/adjoint_processes/replace_elements_and_conditions_for_adjoint_problem_process.h"

//response function
#include "custom_utilities/adjoint_utilities/structural_response_function.h"
#include "custom_utilities/adjoint_utilities/local_stress_response_function.h"
#include "custom_utilities/adjoint_utilities/nodal_displacement_response_function.h"
#include "custom_utilities/adjoint_utilities/strain_energy_response_function.h" 
#include "custom_utilities/adjoint_utilities/eigenfrequency_response_function.h" 



namespace Kratos
{
namespace Python
{

inline
void TransferEigenvector1(
        EigenvectorToSolutionStepVariableTransferUtility& rThisUtil,
        ModelPart& rModelPart,
        int iEigenMode)
{
    rThisUtil.Transfer(rModelPart,iEigenMode);
}

inline
void TransferEigenvector2(
        EigenvectorToSolutionStepVariableTransferUtility& rThisUtil,
        ModelPart& rModelPart,
        int iEigenMode,
        int step)
{
    rThisUtil.Transfer(rModelPart,iEigenMode,step);
}

inline
void CalculateGradient1(
        StructuralResponseFunction& rThisUtil,
        const Condition& rAdjointCondition,
        const Matrix& rAdjointMatrix,
        Vector& rResponseGradient,
        ProcessInfo& rProcessInfo)
{
    rThisUtil.CalculateGradient(rAdjointCondition,rAdjointMatrix,rResponseGradient,rProcessInfo);
}

inline
void CalculateGradient2(
        StructuralResponseFunction& rThisUtil,
        const Element& rAdjointElem,
        const Matrix& rAdjointMatrix,
        Vector& rResponseGradient,
        ProcessInfo& rProcessInfo)
{
    rThisUtil.CalculateGradient(rAdjointElem,rAdjointMatrix,rResponseGradient,rProcessInfo);
}

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

    class_<EigenvectorToSolutionStepVariableTransferUtility>(
                "EigenvectorToSolutionStepVariableTransferUtility")
    .def("Transfer",TransferEigenvector1)
    .def("Transfer",TransferEigenvector2)
    ;

    /// Processes
    class_<ApplyMultipointConstraintsProcess, boost::noncopyable, bases<Process>>("ApplyMultipointConstraintsProcess", init<ModelPart&>())
    .def(init< ModelPart&, Parameters& >())
	.def("AddMasterSlaveRelation", &ApplyMultipointConstraintsProcess::AddMasterSlaveRelationWithNodesAndVariableComponents)
    .def("AddMasterSlaveRelation", &ApplyMultipointConstraintsProcess::AddMasterSlaveRelationWithNodeIdsAndVariableComponents)
	.def("AddMasterSlaveRelation", &ApplyMultipointConstraintsProcess::AddMasterSlaveRelationWithNodesAndVariable)
    .def("AddMasterSlaveRelation", &ApplyMultipointConstraintsProcess::AddMasterSlaveRelationWithNodeIdsAndVariable)
    .def("SetActive", &ApplyMultipointConstraintsProcess::SetActive)      
    .def("PrintData", &ApplyMultipointConstraintsProcess::PrintData);


     class_< OutputPrimalSolutionProcess, bases<Process> >
    ("OutputPrimalSolutionProcess", init<ModelPart&, Parameters&>())
    ;

    class_< InputPrimalSolutionProcess, bases<Process> >
    ("InputPrimalSolutionProcess", init<ModelPart&, Parameters&>())
    ;

    class_<ReplaceElementsAndConditionsForAdjointProblemProcess , bases<Process>, boost::noncopyable >("ReplaceElementsAndConditionsForAdjointProblemProcess",
            init<ModelPart&, Parameters>())
    ;

    //response functions
    class_<StructuralResponseFunction, boost::noncopyable>("StructuralResponseFunction", init<ModelPart&, Parameters&>())
        .def("Initialize", &StructuralResponseFunction::Initialize)
        .def("InitializeSolutionStep", &StructuralResponseFunction::InitializeSolutionStep)
        .def("FinalizeSolutionStep", &StructuralResponseFunction::FinalizeSolutionStep)
        .def("Check", &StructuralResponseFunction::Check)
        .def("Clear", &StructuralResponseFunction::Clear)
        .def("CalculateGradient", CalculateGradient1)
        .def("CalculateGradient", CalculateGradient2)
        .def("CalculateFirstDerivativesGradient",
             &StructuralResponseFunction::CalculateFirstDerivativesGradient)
        .def("CalculateSecondDerivativesGradient",
             &StructuralResponseFunction::CalculateSecondDerivativesGradient)
        .def("CalculateValue", &StructuralResponseFunction::CalculateValue)
        .def("UpdateSensitivities", &StructuralResponseFunction::UpdateSensitivities);  

    class_<LocalStressResponseFunction, bases<StructuralResponseFunction>, boost::noncopyable>
      ("LocalStressResponseFunction", init<ModelPart&, Parameters&>()); 

    class_<NodalDisplacementResponseFunction, bases<StructuralResponseFunction>, boost::noncopyable>
      ("NodalDisplacementResponseFunction", init<ModelPart&, Parameters&>());    

    class_<StrainEnergyResponseFunction, bases<StructuralResponseFunction>, boost::noncopyable>
      ("StrainEnergyResponseFunction", init<ModelPart&, Parameters&>());  

    class_<EigenfrequencyResponseFunction, bases<StructuralResponseFunction>, boost::noncopyable>
      ("EigenfrequencyResponseFunction", init<ModelPart&, Parameters&>()); 

    //for global finite differences
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
        .def("GetStrainEnergy", &FiniteDifferencesUtilities::GetStrainEnergy)
        ;                 

}

}  // namespace Python.  

} // Namespace Kratos

