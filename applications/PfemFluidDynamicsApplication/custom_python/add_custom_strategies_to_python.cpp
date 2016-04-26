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
//#include "custom_strategies/fluid_residual_based_newton_raphson_line_search_implex_strategy.hpp"

#include "custom_strategies/two_step_v_p_strategy.h"
#include "custom_strategies/temporary_predictorcorrector_velocity_bdf_scheme_turbulent.h"

// builder_and_solvers
#include "custom_strategies/custom_builders_and_solvers/residualbased_block_builder_and_solver_periodic.h"

//convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"


//linear solvers
#include "linear_solvers/linear_solver.h"

//schemes
//#include "custom_strategies/custom_schemes/fluid_residual_based_U_wP_static_scheme.hpp"

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
      typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaBaseType;

      //custom types
      // typedef FluidResidualBasedUwPStaticScheme< SparseSpaceType, LocalSpaceType > FluidResidualBasedUwPStaticSchemeType;
      // typedef FluidResidualBasedNewtonRaphsonLineSearchImplexStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > FluidResidualBasedNewtonRaphsonLineSearchImplexStrategyType;

      //********************************************************************
      //*************************STRATEGY CLASSES***************************
      //********************************************************************

    

      //********************************************************************
      //*************************BUILDER AND SOLVER*************************
      //********************************************************************
      


      //********************************************************************
      //*************************SHCHEME CLASSES****************************
      //********************************************************************
    // Static Scheme Type
      // class_< FluidResidualBasedUwPStaticSchemeType,
      // 	      bases< BaseSchemeType >, boost::noncopyable >
      // 	(
      // 	 "FluidResidualBasedUwPStaticScheme", init< >() )
        
      // 	.def("Initialize", &FluidResidualBasedUwPStaticScheme<SparseSpaceType, LocalSpaceType>::Initialize)
      // 	;


      class_< TwoStepVPStrategy< SparseSpaceType,LocalSpaceType, LinearSolverType >, bases<BaseSolvingStrategyType>, boost::noncopyable >
      	("TwoStepVPStrategy",init<ModelPart&,LinearSolverType::Pointer,LinearSolverType::Pointer,bool,bool,double,double,int,int,unsigned int,unsigned int,bool>())
      	//.def(init< ModelPart&, TwoStepVPSolverSettings< SparseSpaceType,LocalSpaceType, LinearSolverType >&, bool >() )
      	//.def(init< ModelPart&, TwoStepVPSolverSettings< SparseSpaceType,LocalSpaceType, LinearSolverType >&, bool, const Kratos::Variable<int>& >() )
      	.def("CalculateReactions",&TwoStepVPStrategy<SparseSpaceType,LocalSpaceType,LinearSolverType>::CalculateReactions)
      	.def("CalculateAccelerations",&TwoStepVPStrategy<SparseSpaceType,LocalSpaceType,LinearSolverType>::CalculateAccelerations)
      	.def("CalculateDisplacements",&TwoStepVPStrategy<SparseSpaceType,LocalSpaceType,LinearSolverType>::CalculateDisplacements)
      	.def("CalculateHistoricalVariables",&TwoStepVPStrategy<SparseSpaceType,LocalSpaceType,LinearSolverType>::CalculateHistoricalVariables)
      	.def("InitializeStressStrain",&TwoStepVPStrategy<SparseSpaceType,LocalSpaceType,LinearSolverType>::InitializeStressStrain)
      	.def("AddIterationStep",&TwoStepVPStrategy<SparseSpaceType,LocalSpaceType,LinearSolverType>::AddIterationStep)
      	.def("ClearExtraIterationSteps",&TwoStepVPStrategy<SparseSpaceType,LocalSpaceType,LinearSolverType>::ClearExtraIterationSteps)
      	;

      //********************************************************************
      //*******************CONVERGENCE CRITERIA CLASSES*********************
      //********************************************************************


     //********************************************************************
    //*************************STRATEGY CLASSES***************************
    //********************************************************************

    // Residual Based Newton-Raphson Line Search Strategy
    // class_< FluidResidualBasedNewtonRaphsonLineSearchImplexStrategyType, 
    //   bases< BaseSolvingStrategyType >, boost::noncopyable >
    //   (
    //    "ResidualBasedNewtonRaphsonLineSearchImplexStrategy",
    //    init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaBaseType::Pointer, int, bool, bool, bool>())
      
    //   .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaBaseType::Pointer, BuilderAndSolverType::Pointer, int, bool, bool, bool >())
    //   .def("SetMaxIterationNumber", &FluidResidualBasedNewtonRaphsonLineSearchImplexStrategyType::SetMaxIterationNumber)
    //   .def("GetMaxIterationNumber", &FluidResidualBasedNewtonRaphsonLineSearchImplexStrategyType::GetMaxIterationNumber)
    //   .def("SetInitializePerformedFlag", &FluidResidualBasedNewtonRaphsonLineSearchImplexStrategyType::SetInitializePerformedFlag)
    //   .def("GetInitializePerformedFlag", &FluidResidualBasedNewtonRaphsonLineSearchImplexStrategyType::GetInitializePerformedFlag)
    //   .def("SetKeepSystemConstantDuringIterations", &FluidResidualBasedNewtonRaphsonLineSearchImplexStrategyType::SetKeepSystemConstantDuringIterations)
    //   .def("GetKeepSystemConstantDuringIterations", &FluidResidualBasedNewtonRaphsonLineSearchImplexStrategyType::GetKeepSystemConstantDuringIterations)
    //   ;
     


    }

  }  // namespace Python.

} // Namespace Kratos

