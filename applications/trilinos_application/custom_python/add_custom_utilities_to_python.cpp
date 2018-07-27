//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//

// System includes

// External includes
#include "Epetra_MpiComm.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"

// Project includes
#include "includes/define.h"
#include "linear_solvers/linear_solver.h"
#include "spaces/ublas_space.h"

// Application includes
#include "trilinos_space.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/trilinos_pointer_wrapper.h"
#include "custom_utilities/trilinos_deactivation_utility.h"
#include "custom_utilities/parallel_fill_communicator.h"
#include "custom_utilities/trilinos_cutting_app.h"
#include "custom_utilities/trilinos_cutting_iso_app.h"
#include "custom_utilities/trilinos_refine_mesh.h"
#include "custom_utilities/trilinos_fractional_step_settings.h"
#include "custom_utilities/trilinos_fractional_step_settings_periodic.h"
#include "custom_utilities/gather_modelpart_utility.h"
#include "custom_utilities/mpi_normal_calculation_utilities.h"
#include "custom_utilities/trilinos_partitioned_fsi_utilities.h"
#include "custom_utilities/trilinos_mvqn_recursive_convergence_accelerator.hpp"

// External includes
#include "../FSIapplication/custom_utilities/convergence_accelerator.hpp"
#include "../FSIapplication/custom_utilities/aitken_convergence_accelerator.hpp"

namespace Kratos
{
namespace Python
{
using namespace pybind11;

typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;
typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
typedef LinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosLinearSolverType;

template <unsigned int TDim>
void AuxiliarUpdateInterfaceValues(
    TrilinosPartitionedFSIUtilities<TrilinosSparseSpaceType, TDim> &dummy,
    ModelPart &rModelPart,
    const Variable<array_1d<double, 3>> &rSolutionVariable,
    AuxiliaryVectorWrapper &rCorrectedGuess)
{
    dummy.UpdateInterfaceValues(rModelPart, rSolutionVariable, rCorrectedGuess.GetReference());
}

template <unsigned int TDim>
void AuxiliarComputeInterfaceResidualVector(
    TrilinosPartitionedFSIUtilities<TrilinosSparseSpaceType, TDim> &dummy,
    ModelPart &rInterfaceModelPart,
    const Variable<array_1d<double, 3>> &rOriginalVariable,
    const Variable<array_1d<double, 3>> &rModifiedVariable,
    AuxiliaryVectorWrapper &rInterfaceResidual)
{
    dummy.ComputeInterfaceResidualVector(rInterfaceModelPart, rOriginalVariable, rModifiedVariable, rInterfaceResidual.GetReference());
}

void AuxiliarUpdateSolution(
    ConvergenceAccelerator<TrilinosSparseSpaceType> &dummy,
    AuxiliaryVectorWrapper &rResidualVector,
    AuxiliaryVectorWrapper &rIterationGuess)
{
    dummy.UpdateSolution(rResidualVector.GetReference(), rIterationGuess.GetReference());
}

void  AddCustomUtilitiesToPython(pybind11::module& m)
{

    class_<TrilinosDeactivationUtility >
    (m,"TrilinosDeactivationUtility")
    .def(init<>() )
    .def("Deactivate", &TrilinosDeactivationUtility::Deactivate )
    .def("Reactivate", &TrilinosDeactivationUtility::Reactivate )
    .def("ReactivateStressFree", &TrilinosDeactivationUtility::ReactivateStressFree )
    .def("ReactivateAll", &TrilinosDeactivationUtility::ReactivateAll )
    .def("Initialize", &TrilinosDeactivationUtility::Initialize )
    ;

    class_<ParallelFillCommunicator >
    (m,"ParallelFillCommunicator")
    .def(init<ModelPart& >() )
    .def("Execute", &ParallelFillCommunicator::Execute )
    .def("PrintDebugInfo", &ParallelFillCommunicator::PrintDebugInfo )
    ;

    class_<TrilinosCuttingApplication>(m,"TrilinosCuttingApplication").def(init< Epetra_MpiComm& >() )
    .def("FindSmallestEdge", &TrilinosCuttingApplication::FindSmallestEdge )
    .def("GenerateCut", &TrilinosCuttingApplication::GenerateCut )
    .def("AddSkinConditions", &TrilinosCuttingApplication::AddSkinConditions )
    .def("AddVariablesToCutModelPart", &TrilinosCuttingApplication::AddVariablesToCutModelPart )
    .def("UpdateCutData", &TrilinosCuttingApplication::UpdateCutData )
    ;

    class_<TrilinosCuttingIsosurfaceApplication >
    (m,"TrilinosCuttingIsosurfaceApplication").def(init< Epetra_MpiComm& >() )
    .def("GenerateScalarVarCut", &TrilinosCuttingIsosurfaceApplication::GenerateVariableCut<double>)
    //.def("GenerateVectorialComponentVarCut", &TrilinosCuttingIsosurfaceApplication::GenerateVectorialComponentVariableCut<VectorComponentAdaptor< array_1d < double, 3 > > >)
    //.def("GenerateVectorialVarCut", &TrilinosCuttingIsosurfaceApplication::GenerateVariableCut< array_1d < double, 3 > >)
    .def("AddSkinConditions", &TrilinosCuttingIsosurfaceApplication::AddSkinConditions)
    .def("UpdateCutData", &TrilinosCuttingIsosurfaceApplication::UpdateCutData)
    .def("DeleteCutData", &TrilinosCuttingIsosurfaceApplication::DeleteCutData)
    ;

    class_<TrilinosRefineMesh>(m,"TrilinosRefineMesh").def(init<ModelPart& , Epetra_MpiComm& >() )
    .def("Local_Refine_Mesh", &TrilinosRefineMesh::Local_Refine_Mesh )
    .def("PrintDebugInfo", &TrilinosRefineMesh::PrintDebugInfo )
    ;

    typedef SolverSettings<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType> BaseSettingsType;
    typedef void (BaseSettingsType::*BuildTurbModelType)(BaseSettingsType::TurbulenceModelLabel const&, TrilinosLinearSolverType::Pointer, const double, const unsigned int);
    typedef void (BaseSettingsType::*PassTurbModelType)(Process::Pointer);
    BuildTurbModelType SetTurbModel_Build = &SolverSettings<TrilinosSparseSpaceType,TrilinosLocalSpaceType,TrilinosLinearSolverType>::SetTurbulenceModel;
    PassTurbModelType SetTurbModel_Pass = &SolverSettings<TrilinosSparseSpaceType,TrilinosLocalSpaceType,TrilinosLinearSolverType>::SetTurbulenceModel;


    class_ < BaseSettingsType >(m,"BaseSettingsType" )
    .def("SetTurbulenceModel",SetTurbModel_Build)
    .def("SetTurbulenceModel",SetTurbModel_Pass)
    ;

    typedef TrilinosFractionalStepSettings<TrilinosSparseSpaceType,TrilinosLocalSpaceType,TrilinosLinearSolverType> TrilinosFSSettingsType;

    enum_<BaseSettingsType::StrategyLabel>(m,"TrilinosStrategyLabel")
    .value("Velocity",BaseSettingsType::Velocity)
    .value("Pressure",BaseSettingsType::Pressure)
    //.value("EddyViscosity",TrilinosFractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::EddyViscosity)
    ;

    enum_<BaseSettingsType::TurbulenceModelLabel>(m,"TrilinosTurbulenceModelLabel")
    .value("SpalartAllmaras",BaseSettingsType::SpalartAllmaras)
    ;

    typedef void (TrilinosFSSettingsType::*SetStrategyByParamsType)(TrilinosFSSettingsType::StrategyLabel const&,TrilinosLinearSolverType::Pointer,const double,const unsigned int);
    SetStrategyByParamsType ThisSetStrategyOverload = &TrilinosFSSettingsType::SetStrategy;

    class_< TrilinosFSSettingsType,BaseSettingsType>(m,"TrilinosFractionalStepSettings")
    .def(init<Epetra_MpiComm&,ModelPart&,unsigned int,unsigned int,bool,bool,bool>())
    .def("SetStrategy",ThisSetStrategyOverload)
    .def("GetStrategy",&TrilinosFSSettingsType::pGetStrategy)
    .def("SetEchoLevel",&TrilinosFSSettingsType::SetEchoLevel)
    ;


    typedef TrilinosFractionalStepSettingsPeriodic<TrilinosSparseSpaceType,TrilinosLocalSpaceType,TrilinosLinearSolverType> TrilinosFSSettingsPeriodicType;

    typedef void (TrilinosFSSettingsPeriodicType::*SetStrategyByParamsPeriodicType)(BaseSettingsType::StrategyLabel const&,TrilinosLinearSolverType::Pointer,const double,const unsigned int);
    SetStrategyByParamsPeriodicType ThatSetStrategyOverload = &TrilinosFSSettingsPeriodicType::SetStrategy;

    class_< TrilinosFSSettingsPeriodicType,BaseSettingsType>
            (m,"TrilinosFractionalStepSettingsPeriodic").def(init<Epetra_MpiComm&,ModelPart&,unsigned int,unsigned int,bool,bool,bool,const Kratos::Variable<int>&>())
    .def("SetStrategy",ThatSetStrategyOverload)
    .def("GetStrategy",&TrilinosFSSettingsPeriodicType::pGetStrategy)
    .def("SetEchoLevel",&TrilinosFSSettingsPeriodicType::SetEchoLevel)
    ;


    class_<GatherModelPartUtility>(m,"GatherModelPartUtility")
        .def(init<int, ModelPart&, int , ModelPart&>() )
      .def("GatherOnMaster",&GatherModelPartUtility::GatherOnMaster<double> )
      .def("GatherOnMaster",&GatherModelPartUtility::GatherOnMaster<array_1d<double,3> > )
      .def("ScatterFromMaster",&GatherModelPartUtility::ScatterFromMaster<double> )
      .def("ScatterFromMaster",&GatherModelPartUtility::ScatterFromMaster<array_1d<double,3> > )
   ;


    class_<MPINormalCalculationUtils, MPINormalCalculationUtils::Pointer>(m,"MPINormalCalculationUtils").def(init<>())
            .def("Check",&MPINormalCalculationUtils::Check)
            .def("OrientFaces",&MPINormalCalculationUtils::OrientFaces)
            .def("CalculateOnSimplex",&MPINormalCalculationUtils::CalculateOnSimplex)
            ;

    typedef PartitionedFSIUtilities<TrilinosSparseSpaceType, 2> BasePartitionedFSIUtilities2DType;
    typedef PartitionedFSIUtilities<TrilinosSparseSpaceType, 3> BasePartitionedFSIUtilities3DType;

    class_<BasePartitionedFSIUtilities2DType, BasePartitionedFSIUtilities2DType::Pointer>(m, "PartitionedFSIUtilities2D");
    class_<BasePartitionedFSIUtilities3DType, BasePartitionedFSIUtilities3DType::Pointer>(m, "PartitionedFSIUtilities3D");

    typedef TrilinosPartitionedFSIUtilities<TrilinosSparseSpaceType,2> TrilinosPartitionedFSIUtilities2DType;
    typedef TrilinosPartitionedFSIUtilities<TrilinosSparseSpaceType,3> TrilinosPartitionedFSIUtilities3DType;

    class_<TrilinosPartitionedFSIUtilities2DType, TrilinosPartitionedFSIUtilities2DType::Pointer, BasePartitionedFSIUtilities2DType>(m, "TrilinosPartitionedFSIUtilities2D")
        .def(init<const Epetra_MpiComm &>())
        .def("GetInterfaceArea", &TrilinosPartitionedFSIUtilities2DType::GetInterfaceArea)
        .def("GetInterfaceResidualSize", &TrilinosPartitionedFSIUtilities2DType::GetInterfaceResidualSize)
        .def("SetUpInterfaceVector", [](TrilinosPartitionedFSIUtilities2DType& self, ModelPart& rModelPart){ 
            return AuxiliaryVectorWrapper(self.SetUpInterfaceVector(rModelPart));})
        .def("UpdateInterfaceValues", &AuxiliarUpdateInterfaceValues<2>)
        .def("ComputeInterfaceResidualVector", &AuxiliarComputeInterfaceResidualVector<2>)
        .def("ComputeFluidInterfaceMeshVelocityResidualNorm", &TrilinosPartitionedFSIUtilities2DType::ComputeFluidInterfaceMeshVelocityResidualNorm)
        .def("ComputeAndPrintFluidInterfaceNorms", &TrilinosPartitionedFSIUtilities2DType::ComputeAndPrintFluidInterfaceNorms)
        .def("ComputeAndPrintStructureInterfaceNorms", &TrilinosPartitionedFSIUtilities2DType::ComputeAndPrintStructureInterfaceNorms)
        .def("CheckCurrentCoordinatesFluid", &TrilinosPartitionedFSIUtilities2DType::CheckCurrentCoordinatesFluid)
        .def("CheckCurrentCoordinatesStructure", &TrilinosPartitionedFSIUtilities2DType::CheckCurrentCoordinatesStructure);

    class_<TrilinosPartitionedFSIUtilities3DType, TrilinosPartitionedFSIUtilities3DType::Pointer, BasePartitionedFSIUtilities3DType>(m, "TrilinosPartitionedFSIUtilities3D")
        .def(init<const Epetra_MpiComm &>())
        .def("GetInterfaceArea", &TrilinosPartitionedFSIUtilities3DType::GetInterfaceArea)
        .def("GetInterfaceResidualSize", &TrilinosPartitionedFSIUtilities3DType::GetInterfaceResidualSize)
        .def("SetUpInterfaceVector", [](TrilinosPartitionedFSIUtilities3DType& self, ModelPart& rModelPart){ 
            return AuxiliaryVectorWrapper(self.SetUpInterfaceVector(rModelPart));})
        .def("UpdateInterfaceValues", &AuxiliarUpdateInterfaceValues<3>)
        .def("ComputeInterfaceResidualVector", &AuxiliarComputeInterfaceResidualVector<3>)
        .def("ComputeFluidInterfaceMeshVelocityResidualNorm", &TrilinosPartitionedFSIUtilities3DType::ComputeFluidInterfaceMeshVelocityResidualNorm)
        .def("ComputeAndPrintFluidInterfaceNorms", &TrilinosPartitionedFSIUtilities3DType::ComputeAndPrintFluidInterfaceNorms)
        .def("ComputeAndPrintStructureInterfaceNorms", &TrilinosPartitionedFSIUtilities3DType::ComputeAndPrintStructureInterfaceNorms)
        .def("CheckCurrentCoordinatesFluid", &TrilinosPartitionedFSIUtilities3DType::CheckCurrentCoordinatesFluid)
        .def("CheckCurrentCoordinatesStructure", &TrilinosPartitionedFSIUtilities3DType::CheckCurrentCoordinatesStructure);

    // Convergence accelerators (from FSIApplication)
    typedef ConvergenceAccelerator<TrilinosSparseSpaceType> TrilinosConvergenceAccelerator;
    typedef AitkenConvergenceAccelerator<TrilinosSparseSpaceType> TrilinosAitkenAccelerator;
    typedef TrilinosMVQNRecursiveJacobianConvergenceAccelerator<TrilinosSparseSpaceType> TrilinosMVQNRecursiveAccelerator;

    // Convergence accelerator base class
    class_< TrilinosConvergenceAccelerator> (m,"TrilinosConvergenceAccelerator").def(init < >())
        .def("Initialize", &TrilinosConvergenceAccelerator::Initialize)
        .def("InitializeSolutionStep", &TrilinosConvergenceAccelerator::InitializeSolutionStep)
        .def("InitializeNonLinearIteration", &TrilinosConvergenceAccelerator::InitializeNonLinearIteration)
        .def("UpdateSolution", AuxiliarUpdateSolution)
        .def("FinalizeNonLinearIteration", &TrilinosConvergenceAccelerator::FinalizeNonLinearIteration)
        .def("FinalizeSolutionStep", &TrilinosConvergenceAccelerator::FinalizeSolutionStep)
        .def("SetEchoLevel", &TrilinosConvergenceAccelerator::SetEchoLevel)
        ;

    class_<TrilinosAitkenAccelerator, TrilinosConvergenceAccelerator>(m,"TrilinosAitkenConvergenceAccelerator")
        .def(init<double>())
        .def(init< Parameters& >())
        ;

    class_< TrilinosMVQNRecursiveAccelerator, TrilinosConvergenceAccelerator>(m,"TrilinosMVQNRecursiveJacobianConvergenceAccelerator")
        .def(init< ModelPart&, const Epetra_MpiComm&, Parameters& >())
        .def(init< ModelPart&, const Epetra_MpiComm&, double, unsigned int >())
        ;

}
}  // namespace Python.

} // Namespace Kratos
