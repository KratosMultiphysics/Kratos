//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:               JMCarbonell $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:               February 2016 $
//   Revision:            $Revision:                     0.0 $
//
//

// System includes 

// External includes 
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/timer.hpp> 


// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_strategies_to_python.h"

#include "spaces/ublas_space.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_strategies/strategies/two_step_v_p_strategy.h"
#include "custom_strategies/strategies/gauss_seidel_linear_strategy.h"
#include "custom_strategies/strategies/explicit_two_step_v_p_strategy.hpp"

//schemes
#include "custom_strategies/schemes/first_order_forward_euler_scheme.hpp" 

// builder_and_solvers

//convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"


//linear solvers
#include "linear_solvers/linear_solver.h"

//schemes

namespace Kratos
{

  namespace Python
  {		
    using namespace boost::python;

    void  AddCustomStrategiesToPython()
    {
      typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
      typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

      //base types
      typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
      typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
      typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;
      typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;
      //typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaBaseType;
      typedef FirstOrderForwardEulerScheme< SparseSpaceType, LocalSpaceType >  FirstOrderForwardEulerSchemeType;

      //custom strategy types
      typedef ExplicitTwoStepVPStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ExplicitStrategyType;


      //********************************************************************
      //*************************SHCHEME CLASSES****************************
      //********************************************************************

      class_< TwoStepVPStrategy< SparseSpaceType,LocalSpaceType, LinearSolverType >, bases<BaseSolvingStrategyType>, boost::noncopyable >
      	("TwoStepVPStrategy",init<ModelPart&,LinearSolverType::Pointer,LinearSolverType::Pointer,bool,double,double,int,unsigned int,unsigned int>())
      	.def("CalculateAccelerations",&TwoStepVPStrategy<SparseSpaceType,LocalSpaceType,LinearSolverType>::CalculateAccelerations)
      	.def("CalculateDisplacements",&TwoStepVPStrategy<SparseSpaceType,LocalSpaceType,LinearSolverType>::CalculateDisplacements)
      	// .def("InitializeStressStrain",&TwoStepVPStrategy<SparseSpaceType,LocalSpaceType,LinearSolverType>::InitializeStressStrain)
	;

      class_< GaussSeidelLinearStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >, bases< BaseSolvingStrategyType >, boost::noncopyable >
	("GaussSeidelLinearStrategy",init < ModelPart& ,  BaseSchemeType::Pointer, LinearSolverType::Pointer, BuilderAndSolverType::Pointer, bool, bool  >())
	.def("GetResidualNorm", &GaussSeidelLinearStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::GetResidualNorm)
	.def("SetBuilderAndSolver", &GaussSeidelLinearStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SetBuilderAndSolver)
	;

      class_< ExplicitStrategyType, 
	      bases< BaseSolvingStrategyType >, boost::noncopyable >
	(
	 "ExplicitStrategy",
	 init < ModelPart&, BaseSchemeType::Pointer,  LinearSolverType::Pointer, bool, bool, bool >())
      
	.def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer,  bool, bool, bool >())
	.def("SetInitializePerformedFlag", &ExplicitStrategyType::SetInitializePerformedFlag)
	.def("GetInitializePerformedFlag", &ExplicitStrategyType::GetInitializePerformedFlag)
	;

      // Explicit scheme: Central differences 
      class_< FirstOrderForwardEulerSchemeType,
	      bases< BaseSchemeType >,  boost::noncopyable >
	(
	 "FirstOrderForwardEulerScheme", init< const double, const double, const double, const bool >() )

	.def("Initialize", &FirstOrderForwardEulerScheme<SparseSpaceType, LocalSpaceType>::Initialize)
	;
    }

  }  // namespace Python.

} // Namespace Kratos

