//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
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
#include "custom_strategies/strategies/residual_based_newton_raphson_line_search_implex_strategy.hpp"

//convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"


//linear solvers
#include "linear_solvers/linear_solver.h"

//schemes
#include "custom_strategies/schemes/residual_based_static_scheme.hpp"
#include "custom_strategies/schemes/residual_based_newmark_scheme.hpp"
#include "custom_strategies/schemes/residual_based_bossak_scheme.hpp"
#include "custom_strategies/schemes/residual_based_contact_bossak_scheme.hpp"
#include "custom_strategies/schemes/residual_based_U_wP_static_scheme.hpp"

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

      //custom strategy types
      typedef ResidualBasedUwPStaticScheme< SparseSpaceType, LocalSpaceType > ResidualBasedUwPStaticSchemeType;
      typedef ResidualBasedNewtonRaphsonLineSearchImplexStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedNewtonRaphsonLineSearchImplexStrategyType;

      //custom scheme types
      typedef ResidualBasedStaticScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedStaticSchemeType;
      typedef ResidualBasedNewmarkScheme< SparseSpaceType, LocalSpaceType > ResidualBasedNewmarkSchemeType;
      typedef ResidualBasedBossakScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedBossakSchemeType;
      typedef ResidualBasedContactBossakScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedContactBossakSchemeType;    

    

      //********************************************************************
      //*************************BUILDER AND SOLVER*************************
      //********************************************************************
      


      //********************************************************************
      //*************************SHCHEME CLASSES****************************
      //********************************************************************

      // Static Scheme Type
      class_< ResidualBasedStaticSchemeType,
	      bases< BaseSchemeType >, boost::noncopyable >
	(
	 "ResidualBasedStaticScheme", init< >() )
	
	.def("Initialize", &ResidualBasedStaticScheme<SparseSpaceType, LocalSpaceType>::Initialize)
	;
      
      // Residual Based Newmark Scheme Type
      class_< ResidualBasedNewmarkSchemeType,
	      bases< BaseSchemeType >, boost::noncopyable >
	(
	 "ResidualBasedNewmarkScheme", init< double >() )
	
	.def("Initialize", &ResidualBasedNewmarkScheme<SparseSpaceType, LocalSpaceType>::Initialize)
	;
      
      // Residual Based Bossak Scheme Type
      class_< ResidualBasedBossakSchemeType,
	      bases< BaseSchemeType >,  boost::noncopyable >
	(
	 "ResidualBasedBossakScheme", init< double , double >() )
	
	.def("Initialize", &ResidualBasedBossakScheme<SparseSpaceType, LocalSpaceType>::Initialize)
	;
      

      // Residual Based Bossak Scheme Type
      class_< ResidualBasedContactBossakSchemeType,
	      bases< BaseSchemeType >,  boost::noncopyable >
	(
	 "ResidualBasedContactBossakScheme", init< double , double >() )
	
	.def("Initialize", &ResidualBasedContactBossakScheme<SparseSpaceType, LocalSpaceType>::Initialize)
	;
      
      // Static Scheme Type
      class_< ResidualBasedUwPStaticSchemeType,
	      bases< BaseSchemeType >, boost::noncopyable >
	(
	 "ResidualBasedUwPStaticScheme", init< >() )
        
	.def("Initialize", &ResidualBasedUwPStaticScheme<SparseSpaceType, LocalSpaceType>::Initialize)
	;
      
      
      //********************************************************************
      //*******************CONVERGENCE CRITERIA CLASSES*********************
      //********************************************************************
      
      
      //********************************************************************
      //*************************STRATEGY CLASSES***************************
      //********************************************************************
      
      // Residual Based Newton-Raphson Line Search Strategy
      class_< ResidualBasedNewtonRaphsonLineSearchImplexStrategyType, 
	      bases< BaseSolvingStrategyType >, boost::noncopyable >
	(
	 "ResidualBasedNewtonRaphsonLineSearchImplexStrategy",
	 init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaBaseType::Pointer, int, bool, bool, bool>())
	
	.def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaBaseType::Pointer, BuilderAndSolverType::Pointer, int, bool, bool, bool >())
	.def("SetMaxIterationNumber", &ResidualBasedNewtonRaphsonLineSearchImplexStrategyType::SetMaxIterationNumber)
	.def("GetMaxIterationNumber", &ResidualBasedNewtonRaphsonLineSearchImplexStrategyType::GetMaxIterationNumber)
	.def("SetInitializePerformedFlag", &ResidualBasedNewtonRaphsonLineSearchImplexStrategyType::SetInitializePerformedFlag)
	.def("GetInitializePerformedFlag", &ResidualBasedNewtonRaphsonLineSearchImplexStrategyType::GetInitializePerformedFlag)
	.def("SetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonLineSearchImplexStrategyType::SetKeepSystemConstantDuringIterations)
	.def("GetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonLineSearchImplexStrategyType::GetKeepSystemConstantDuringIterations)
	;
                             
    }

  }  // namespace Python.

} // Namespace Kratos

