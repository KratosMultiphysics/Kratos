//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:               JMCarbonell $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:               February 2016 $
//   Revision:            $Revision:                     0.0 $
//
//

// System includes 
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// External includes 

// Project includes
#include "includes/node.h"
#include "includes/define.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "utilities/openmp_utils.h"

//Application includes
#include "custom_python/add_custom_utilities_to_python.h"
//#include "custom_strategies/custom_builders_and_solvers/residual_based_builder_and_solver.hpp"

//Utilities
// #include "custom_utilities/boundary_normals_calculation_utilities.hpp"
// #include "custom_utilities/modeler_utilities.hpp"
// #include "custom_utilities/contact_domain_utilities.hpp"

#include "custom_utilities/two_step_v_p_settings.h"
#include "custom_utilities/two_step_v_p_settings_periodic.h"

namespace Kratos
{
	
  namespace Python
  {
    
    void  AddCustomUtilitiesToPython()
    {

      using namespace boost::python;

      
      typedef UblasSpace<double, CompressedMatrix, Vector>    SparseSpaceType;
      typedef UblasSpace<double, Matrix, Vector>               LocalSpaceType;
      typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

      typedef Scheme< SparseSpaceType, LocalSpaceType >            SchemeType;
      typedef SchemeType::Pointer                           SchemePointerType;

      typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType >                 BuilderAndSolverType;
      typedef BuilderAndSolverType::Pointer        BuilderAndSolverPointerType;

     //***************DOMAIN SET**************//
      // class_< ModelerUtilities, boost::noncopyable > ("ModelerUtilities", init<>())
      // 	.def("SetDomainLabels",&ModelerUtilities::SetDomainLabels)
      // 	;


      //***************NORMALS**************//

      // // This is required to recognize the different overloads 
      // typedef  void (BoundaryNormalsCalculationUtilities::*CalculateMeshBoundaryNormals)(ModelPart&, int, int);
      // typedef  void (BoundaryNormalsCalculationUtilities::*CalculateMeshUnitBoundaryNormals)(ModelPart&, int, int); 


      // enum_<TwoStepVPSolverSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::TurbulenceModelLabel>("TwoStepVPTurbulenceModelLabel")
      // 	.value("SpalartAllmaras",TwoStepVPSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::SpalartAllmaras)
      // 	;

      // typedef TwoStepVPSolverSettings<SparseSpaceType,LocalSpaceType,LinearSolverType> BaseSettingsType;
   
      // class_ < BaseSettingsType, boost::noncopyable >
      // 	( "TwoStepVPBaseSettingsType",no_init );

      // enum_<TwoStepVPSolverSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::StrategyLabel>("TwoStepVPStrategyLabel")
      // 	.value("Velocity",TwoStepVPSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::Velocity)
      // 	.value("Pressure",TwoStepVPSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::Pressure)
      // 	//.value("EddyViscosity",TwoStepVPSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::EddyViscosity)
      // 	;

      // typedef void (TwoStepVPSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::*SetStrategyByParamsType)(TwoStepVPSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::StrategyLabel const&,LinearSolverType::Pointer,const double,const unsigned int);
      // typedef void (TwoStepVPSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::*BuildTurbModelType)(TwoStepVPSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::TurbulenceModelLabel const&, LinearSolverType::Pointer, const double, const unsigned int);
      // typedef void (TwoStepVPSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::*PassTurbModelType)(Process::Pointer);
      // SetStrategyByParamsType ThisSetStrategyOverload = &TwoStepVPSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetStrategy;
      // BuildTurbModelType SetTurbModel_Build = &TwoStepVPSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetTurbulenceModel;
      // PassTurbModelType SetTurbModel_Pass = &TwoStepVPSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetTurbulenceModel;

      // class_< TwoStepVPSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>,bases<BaseSettingsType>, boost::noncopyable>
      // 	("TwoStepVPSettings",init<ModelPart&,unsigned int,unsigned int,bool,bool,bool>())
      // 	.def("SetStrategy",ThisSetStrategyOverload)
      // 	.def("SetTurbulenceModel",SetTurbModel_Build)
      // 	.def("SetTurbulenceModel",SetTurbModel_Pass)
      // 	.def("GetStrategy",&TwoStepVPSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::pGetStrategy)
      // 	.def("SetEchoLevel",&TwoStepVPSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetEchoLevel)
      // 	;
    
      // class_< TwoStepVPSettingsPeriodic<SparseSpaceType,LocalSpaceType,LinearSolverType>,bases<BaseSettingsType>, boost::noncopyable>
      // 	("TwoStepVPSettingsPeriodic",init<ModelPart&,unsigned int,unsigned int,bool,bool,bool,const Kratos::Variable<int>&>())
      // 	.def("SetStrategy",&TwoStepVPSettingsPeriodic<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetStrategy)
      // 	.def("GetStrategy",&TwoStepVPSettingsPeriodic<SparseSpaceType,LocalSpaceType,LinearSolverType>::pGetStrategy)
      // 	.def("SetEchoLevel",&TwoStepVPSettingsPeriodic<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetEchoLevel)
      // 	;




//      CalculateMeshBoundaryNormals          CalculateMeshNormals     = &BoundaryNormalsCalculationUtilities::CalculateMeshBoundaryNormals;
//      CalculateMeshUnitBoundaryNormals      CalculateMeshUnitNormals = &BoundaryNormalsCalculationUtilities::CalculateMeshUnitBoundaryNormals;
      
//      class_<BoundaryNormalsCalculationUtilities > ("BoundaryNormalsCalculation", init<>())
//	.def("CalculateBoundaryNormals", &BoundaryNormalsCalculationUtilities::CalculateBoundaryNormals)
//	.def("CalculateBoundaryUnitNormals", &BoundaryNormalsCalculationUtilities::CalculateUnitBoundaryNormals)
//	.def("CalculateMeshBoundaryNormals", CalculateMeshNormals)
//	.def("CalculateMeshBoundaryUnitNormals", CalculateMeshUnitNormals)
//	;
    }

  }  // namespace Python.

} // Namespace Kratos

