/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last modified by:    $Author: janosch $
//   Date:                $Date: 2008-04-28 16:19:49 $
//   Revision:            $Revision: 1.2 $
//
//


// System includes

// External includes
#include <boost/python.hpp>
#include "Epetra_MpiComm.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"

// Project includes
#include "includes/define.h"
#include "linear_solvers/linear_solver.h"
#include "spaces/ublas_space.h"

// Application includes
#include "trilinos_space.h"
#include "custom_utilities/trilinos_deactivation_utility.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/parallel_fill_communicator.h"
#include "custom_utilities/trilinos_cutting_app.h"
#include "custom_utilities/trilinos_cutting_iso_app.h"
#include "custom_utilities/trilinos_refine_mesh.h"
#include "custom_utilities/trilinos_fractional_step_settings.h"
#include "custom_utilities/trilinos_fractional_step_settings_periodic.h"
#include "custom_utilities/gather_modelpart_utility.h"
#include "custom_utilities/mpi_normal_calculation_utilities.h"

namespace Kratos
{
namespace Python
{
using namespace boost::python;

void  AddCustomUtilitiesToPython()
{
    typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;
    typedef LinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosLinearSolverType;

    class_<TrilinosDeactivationUtility, boost::noncopyable >
    ("TrilinosDeactivationUtility",
     init<>() )
    .def("Deactivate", &TrilinosDeactivationUtility::Deactivate )
    .def("Reactivate", &TrilinosDeactivationUtility::Reactivate )
    .def("ReactivateStressFree", &TrilinosDeactivationUtility::ReactivateStressFree )
    .def("ReactivateAll", &TrilinosDeactivationUtility::ReactivateAll )
    .def("Initialize", &TrilinosDeactivationUtility::Initialize )
    ;

    class_<ParallelFillCommunicator, boost::noncopyable >
    ("ParallelFillCommunicator",
     init<ModelPart& >() )
    .def("Execute", &ParallelFillCommunicator::Execute )
    .def("PrintDebugInfo", &ParallelFillCommunicator::PrintDebugInfo )
    ;

    class_<TrilinosCuttingApplication, boost::noncopyable >
    ("TrilinosCuttingApplication",
     init< Epetra_MpiComm& >() )
    .def("FindSmallestEdge", &TrilinosCuttingApplication::FindSmallestEdge )
    .def("GenerateCut", &TrilinosCuttingApplication::GenerateCut )
    .def("AddSkinConditions", &TrilinosCuttingApplication::AddSkinConditions )
    .def("UpdateCutData", &TrilinosCuttingApplication::UpdateCutData )
    ;

    class_<TrilinosCuttingIsosurfaceApplication,boost::noncopyable >
    ("TrilinosCuttingIsosurfaceApplication",
     init< Epetra_MpiComm& >() )
    .def("GenerateScalarVarCut", &TrilinosCuttingIsosurfaceApplication::GenerateVariableCut<double>)
    //.def("GenerateVectorialComponentVarCut", &TrilinosCuttingIsosurfaceApplication::GenerateVectorialComponentVariableCut<VectorComponentAdaptor< array_1d < double, 3 > > >)
    //.def("GenerateVectorialVarCut", &TrilinosCuttingIsosurfaceApplication::GenerateVariableCut< array_1d < double, 3 > >)
    .def("AddSkinConditions", &TrilinosCuttingIsosurfaceApplication::AddSkinConditions)
    .def("UpdateCutData", &TrilinosCuttingIsosurfaceApplication::UpdateCutData)
    .def("DeleteCutData", &TrilinosCuttingIsosurfaceApplication::DeleteCutData)
    ;

    class_<TrilinosRefineMesh, boost::noncopyable >
    ("TrilinosRefineMesh",
     init<ModelPart& , Epetra_MpiComm& >() )
    .def("Local_Refine_Mesh", &TrilinosRefineMesh::Local_Refine_Mesh )
    .def("PrintDebugInfo", &TrilinosRefineMesh::PrintDebugInfo )
    ;

    typedef SolverSettings<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType> BaseSettingsType;
    typedef void (BaseSettingsType::*BuildTurbModelType)(BaseSettingsType::TurbulenceModelLabel const&, TrilinosLinearSolverType::Pointer, const double, const unsigned int);
    typedef void (BaseSettingsType::*PassTurbModelType)(Process::Pointer);
    BuildTurbModelType SetTurbModel_Build = &SolverSettings<TrilinosSparseSpaceType,TrilinosLocalSpaceType,TrilinosLinearSolverType>::SetTurbulenceModel;
    PassTurbModelType SetTurbModel_Pass = &SolverSettings<TrilinosSparseSpaceType,TrilinosLocalSpaceType,TrilinosLinearSolverType>::SetTurbulenceModel;


    class_ < BaseSettingsType, boost::noncopyable >
    ( "BaseSettingsType",no_init )
    .def("SetTurbulenceModel",SetTurbModel_Build)
    .def("SetTurbulenceModel",SetTurbModel_Pass)
    ;

    typedef TrilinosFractionalStepSettings<TrilinosSparseSpaceType,TrilinosLocalSpaceType,TrilinosLinearSolverType> TrilinosFSSettingsType;

    enum_<BaseSettingsType::StrategyLabel>("TrilinosStrategyLabel")
    .value("Velocity",BaseSettingsType::Velocity)
    .value("Pressure",BaseSettingsType::Pressure)
    //.value("EddyViscosity",TrilinosFractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::EddyViscosity)
    ;

    enum_<BaseSettingsType::TurbulenceModelLabel>("TrilinosTurbulenceModelLabel")
    .value("SpalartAllmaras",BaseSettingsType::SpalartAllmaras)
    ;
    
    typedef void (TrilinosFSSettingsType::*SetStrategyByParamsType)(TrilinosFSSettingsType::StrategyLabel const&,TrilinosLinearSolverType::Pointer,const double,const unsigned int);
    SetStrategyByParamsType ThisSetStrategyOverload = &TrilinosFSSettingsType::SetStrategy;

    class_< TrilinosFSSettingsType,bases<BaseSettingsType>, boost::noncopyable>
            ("TrilinosFractionalStepSettings",init<Epetra_MpiComm&,ModelPart&,unsigned int,unsigned int,bool,bool,bool>())
    .def("SetStrategy",ThisSetStrategyOverload)
    .def("GetStrategy",&TrilinosFSSettingsType::pGetStrategy)
    .def("SetEchoLevel",&TrilinosFSSettingsType::SetEchoLevel)
    ;


    typedef TrilinosFractionalStepSettingsPeriodic<TrilinosSparseSpaceType,TrilinosLocalSpaceType,TrilinosLinearSolverType> TrilinosFSSettingsPeriodicType;

    typedef void (TrilinosFSSettingsPeriodicType::*SetStrategyByParamsPeriodicType)(BaseSettingsType::StrategyLabel const&,TrilinosLinearSolverType::Pointer,const double,const unsigned int);
    SetStrategyByParamsPeriodicType ThatSetStrategyOverload = &TrilinosFSSettingsPeriodicType::SetStrategy;

    class_< TrilinosFSSettingsPeriodicType,bases<BaseSettingsType>, boost::noncopyable>
            ("TrilinosFractionalStepSettingsPeriodic",init<Epetra_MpiComm&,ModelPart&,unsigned int,unsigned int,bool,bool,bool,const Kratos::Variable<int>&>())
    .def("SetStrategy",ThatSetStrategyOverload)
    .def("GetStrategy",&TrilinosFSSettingsPeriodicType::pGetStrategy)
    .def("SetEchoLevel",&TrilinosFSSettingsPeriodicType::SetEchoLevel)
    ;
    
    
    class_<GatherModelPartUtility, boost::noncopyable >
    ("GatherModelPartUtility",
     init<int, ModelPart&, int , ModelPart&>() )
      .def("GatherOnMaster",&GatherModelPartUtility::GatherOnMaster<double> )
      .def("GatherOnMaster",&GatherModelPartUtility::GatherOnMaster<array_1d<double,3> > )
      .def("ScatterFromMaster",&GatherModelPartUtility::ScatterFromMaster<double> )
      .def("ScatterFromMaster",&GatherModelPartUtility::ScatterFromMaster<array_1d<double,3> > )
   ;


    class_<MPINormalCalculationUtils, MPINormalCalculationUtils::Pointer, boost::noncopyable > ("MPINormalCalculationUtils",init<>())
            .def("Check",&MPINormalCalculationUtils::Check)
            .def("OrientFaces",&MPINormalCalculationUtils::OrientFaces)
            .def("CalculateOnSimplex",&MPINormalCalculationUtils::CalculateOnSimplex)
            ;


}
}  // namespace Python.

} // Namespace Kratos

