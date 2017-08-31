//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/timer.hpp>

// Project includes
#include "includes/define.h"
#include "containers/flags.h"
#include "custom_python/add_custom_strategies_to_python.h"

#include "spaces/ublas_space.h"

//strategies
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

#include "custom_strategies/strategies/component_wise_newton_raphson_strategy.hpp"
#include "custom_strategies/strategies/residual_based_newton_raphson_line_search_strategy.hpp"
#include "custom_strategies/strategies/residual_based_newton_raphson_line_search_implex_strategy.hpp"
#include "custom_strategies/strategies/explicit_strategy.hpp" 
#include "custom_strategies/strategies/eigensolver_strategy.hpp" 
#include "custom_strategies/strategies/explicit_hamilton_strategy.hpp"

//builders and solvers
#include "custom_strategies/builders_and_solvers/component_wise_builder_and_solver.hpp"
#include "custom_strategies/builders_and_solvers/explicit_hamilton_builder_and_solver.hpp"

//convergence criteria
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "custom_strategies/convergence_criteria/displacement_convergence_criterion.hpp"
#include "custom_strategies/convergence_criteria/component_wise_residual_convergence_criterion.hpp"

//schemes
#include "custom_strategies/schemes/component_wise_bossak_scheme.hpp"
#include "custom_strategies/schemes/eigensolver_dynamic_scheme.hpp" 
 
#include "custom_strategies/schemes/explicit_central_differences_scheme.hpp" 
#include "custom_strategies/schemes/residual_based_rotation_newmark_scheme.hpp"
#include "custom_strategies/schemes/residual_based_rotation_simo_scheme.hpp"
#include "custom_strategies/schemes/residual_based_rotation_emc_scheme.hpp"
#include "custom_strategies/schemes/explicit_hamilton_scheme.hpp"

//linear solvers
#include "linear_solvers/linear_solver.h"



namespace Kratos
{

  namespace Python
  {
    using namespace boost::python;

    void  AddCustomStrategiesToPython()
    {
      //base types
      typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
      typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
      typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
      typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
      typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;
      typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;
      typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaType;
      typedef BuilderAndSolverType::Pointer BuilderAndSolverPointer;


      //custom strategy types
      typedef ResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedNewtonRaphsonStrategyType;
      typedef ComponentWiseNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ComponentWiseNewtonRaphsonStrategyType;
      typedef ResidualBasedNewtonRaphsonLineSearchStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedNewtonRaphsonLineSearchStrategyType;
      typedef ResidualBasedNewtonRaphsonLineSearchImplexStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedNewtonRaphsonLineSearchImplexStrategyType;
      typedef ExplicitStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ExplicitStrategyType;
      typedef ExplicitHamiltonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ExplicitHamiltonStrategyType;
      typedef EigensolverStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > EigensolverStrategyType;


      //custom builder_and_solver types
      typedef ComponentWiseBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > ComponentWiseBuilderAndSolverType;
      typedef ExplicitHamiltonBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > ExplicitHamiltonBuilderAndSolverType;

      //custom scheme types
      typedef ComponentWiseBossakScheme< SparseSpaceType, LocalSpaceType >  ComponentWiseBossakSchemeType;     
      typedef ExplicitCentralDifferencesScheme< SparseSpaceType, LocalSpaceType >  ExplicitCentralDifferencesSchemeType;
      typedef ResidualBasedRotationNewmarkScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedRotationNewmarkSchemeType;
      typedef ResidualBasedRotationSimoScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedRotationSimoSchemeType;
      typedef ResidualBasedRotationEMCScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedRotationEMCSchemeType;
      typedef ExplicitHamiltonScheme< SparseSpaceType, LocalSpaceType >  ExplicitHamiltonSchemeType;
      typedef EigensolverDynamicScheme< SparseSpaceType, LocalSpaceType > EigensolverDynamicSchemeType;

      //custom convergence criterion types
      typedef DisplacementConvergenceCriterion< SparseSpaceType,  LocalSpaceType > DisplacementConvergenceCriterionType;
      typedef ComponentWiseResidualConvergenceCriterion< SparseSpaceType,  LocalSpaceType > ComponentWiseResidualConvergenceCriterionType;

      //********************************************************************
      //*************************STRATEGY CLASSES***************************
      //********************************************************************

      // Solid Mechanics Explicit Strategy
      class_< ExplicitStrategyType, 
	      bases< BaseSolvingStrategyType >, boost::noncopyable >
	(
	 "ExplicitStrategy",
	 init < ModelPart&, BaseSchemeType::Pointer,  LinearSolverType::Pointer, bool, bool, bool >())
      
	.def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer,  bool, bool, bool >())
	.def("SetInitializePerformedFlag", &ExplicitStrategyType::SetInitializePerformedFlag)
	.def("GetInitializePerformedFlag", &ExplicitStrategyType::GetInitializePerformedFlag)
	;

        
      // Component Wise Newton-Raphson Strategy
      class_< ComponentWiseNewtonRaphsonStrategyType, 
	      bases< ResidualBasedNewtonRaphsonStrategyType >, boost::noncopyable >
	(
	 "ComponentWiseNewtonRaphsonStrategy",
	 init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, int, bool, bool, bool >())
      
	.def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, int, bool, bool, bool >())
	.def("SetMaxIterationNumber", &ComponentWiseNewtonRaphsonStrategyType::SetMaxIterationNumber)
	.def("GetMaxIterationNumber", &ComponentWiseNewtonRaphsonStrategyType::GetMaxIterationNumber)
	.def("SetKeepSystemConstantDuringIterations", &ComponentWiseNewtonRaphsonStrategyType::SetKeepSystemConstantDuringIterations)
	.def("GetKeepSystemConstantDuringIterations", &ComponentWiseNewtonRaphsonStrategyType::GetKeepSystemConstantDuringIterations)

	;
  

      // Residual Based Newton-Raphson Line Search Strategy
      class_< ResidualBasedNewtonRaphsonLineSearchStrategyType, 
	      bases< BaseSolvingStrategyType >, boost::noncopyable >
	(
	 "ResidualBasedNewtonRaphsonLineSearchStrategy",
	 init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, int, bool, bool, bool>())
      
	.def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, int, bool, bool, bool >())
	.def("SetMaxIterationNumber", &ResidualBasedNewtonRaphsonLineSearchStrategyType::SetMaxIterationNumber)
	.def("GetMaxIterationNumber", &ResidualBasedNewtonRaphsonLineSearchStrategyType::GetMaxIterationNumber)
	.def("SetInitializePerformedFlag", &ResidualBasedNewtonRaphsonLineSearchStrategyType::SetInitializePerformedFlag)
	.def("GetInitializePerformedFlag", &ResidualBasedNewtonRaphsonLineSearchStrategyType::GetInitializePerformedFlag)
	.def("SetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonLineSearchStrategyType::SetKeepSystemConstantDuringIterations)
	.def("GetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonLineSearchStrategyType::GetKeepSystemConstantDuringIterations)
	.def("SetFinalizeSolutionStepFlag", &ResidualBasedNewtonRaphsonLineSearchStrategyType::SetFinalizeSolutionStepFlag)
	.def("GetFinalizeSolutionStepFlag", &ResidualBasedNewtonRaphsonLineSearchStrategyType::GetFinalizeSolutionStepFlag)
	;
     
      
      // Residual Based Newton-Raphson Line Search Implex Strategy
      class_< ResidualBasedNewtonRaphsonLineSearchImplexStrategyType, 
	      bases< BaseSolvingStrategyType >, boost::noncopyable >
	(
	 "ResidualBasedNewtonRaphsonLineSearchImplexStrategy",
	 init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, int, bool, bool, bool>())
	
	.def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, int, bool, bool, bool >())
	.def("SetMaxIterationNumber", &ResidualBasedNewtonRaphsonLineSearchImplexStrategyType::SetMaxIterationNumber)
	.def("GetMaxIterationNumber", &ResidualBasedNewtonRaphsonLineSearchImplexStrategyType::GetMaxIterationNumber)
	.def("SetInitializePerformedFlag", &ResidualBasedNewtonRaphsonLineSearchImplexStrategyType::SetInitializePerformedFlag)
	.def("GetInitializePerformedFlag", &ResidualBasedNewtonRaphsonLineSearchImplexStrategyType::GetInitializePerformedFlag)
	.def("SetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonLineSearchImplexStrategyType::SetKeepSystemConstantDuringIterations)
	.def("GetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonLineSearchImplexStrategyType::GetKeepSystemConstantDuringIterations)
	;

      // Explicit Hamilton Estrategy for Explicit Beam solution
      class_< ExplicitHamiltonStrategyType, 
	      bases< BaseSolvingStrategyType >, boost::noncopyable >
	(
	 "ExplicitHamiltonStrategy",
	 init < ModelPart&, BaseSchemeType::Pointer,  LinearSolverType::Pointer, bool, bool, bool>())
	
	.def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer,  bool, bool, bool >())
	.def("SetInitializePerformedFlag", &ExplicitHamiltonStrategyType::SetInitializePerformedFlag)
	.def("GetInitializePerformedFlag", &ExplicitHamiltonStrategyType::GetInitializePerformedFlag)
	;

      // Eigensolver Scheme Type
      class_< EigensolverDynamicSchemeType, EigensolverDynamicSchemeType::Pointer,
	      bases< BaseSchemeType >, boost::noncopyable >
	(
	 "EigensolverDynamicScheme",
	 init<>())
	;

      
      //********************************************************************
      //*******************BUILDER AND SOLVER CLASSES***********************
      //********************************************************************


      // Component Wise Builder and Solver
      class_< ComponentWiseBuilderAndSolverType, bases<BuilderAndSolverType>, boost::noncopyable > 
	(
	 "ComponentWiseBuilderAndSolver", init< LinearSolverType::Pointer > ()
	 );


      // Residual Based Builder and Solver
      class_< ExplicitHamiltonBuilderAndSolverType, bases<BuilderAndSolverType>, boost::noncopyable > 
	(
	 "ExplicitHamiltonBuilderAndSolver", init< LinearSolverType::Pointer > ()
	 );

      
      //********************************************************************
      //*************************SHCHEME CLASSES****************************
      //********************************************************************


      // Component Wise Bossak Scheme Type
      class_< ComponentWiseBossakSchemeType,
	      bases< BaseSchemeType >,  boost::noncopyable >
	(
	 "ComponentWiseBossakScheme", init< double >() )

	.def("Initialize", &ComponentWiseBossakScheme<SparseSpaceType, LocalSpaceType>::Initialize)
	;


      // Explicit scheme: Central differences 
      class_< ExplicitCentralDifferencesSchemeType,
	      bases< BaseSchemeType >,  boost::noncopyable >
	(
	 "ExplicitCentralDifferencesScheme", init< const double, const double, const double, const bool >() )

	.def("Initialize", &ExplicitCentralDifferencesScheme<SparseSpaceType, LocalSpaceType>::Initialize)
	;

      // Residual Based Rotational Newmark Scheme Type
      class_< ResidualBasedRotationNewmarkSchemeType,
	      bases< BaseSchemeType >,  boost::noncopyable >
	(
	 "ResidualBasedRotationNewmarkScheme", init< double, double >() )
	
	.def("Initialize", &ResidualBasedRotationNewmarkScheme<SparseSpaceType, LocalSpaceType>::Initialize)
	;

      // Residual Based Rotational Simo Scheme Type
      class_< ResidualBasedRotationSimoSchemeType,
	      bases< BaseSchemeType >,  boost::noncopyable >
	(
	 "ResidualBasedRotationSimoScheme", init< double, double >() )
	
	.def("Initialize", &ResidualBasedRotationSimoScheme<SparseSpaceType, LocalSpaceType>::Initialize)
	;

      // Residual Based Rotational EMC Scheme Type
      class_< ResidualBasedRotationEMCSchemeType,
	      bases< BaseSchemeType >,  boost::noncopyable >
	(
	 "ResidualBasedRotationEMCScheme", init< double >() )
	
	.def("Initialize", &ResidualBasedRotationEMCScheme<SparseSpaceType, LocalSpaceType>::Initialize)
	;

      // Explicit Hamilton Scheme Type
      class_< ExplicitHamiltonSchemeType,
	      bases< BaseSchemeType >,  boost::noncopyable >
	(
	 "ExplicitHamiltonScheme", init< double, double, double, bool >() )
	
	.def("Initialize", &ExplicitHamiltonScheme<SparseSpaceType, LocalSpaceType>::Initialize)
	;
      
      //********************************************************************
      //*******************CONVERGENCE CRITERIA CLASSES*********************
      //********************************************************************

      // Displacement Convergence Criterion
      class_< DisplacementConvergenceCriterionType,
	      bases< ConvergenceCriteriaType >, boost::noncopyable >
	(
	 "DisplacementConvergenceCriterion", 
	 init<double, double >())
	.def(init<double, double >())
	.def("SetEchoLevel", &DisplacementConvergenceCriterionType::SetEchoLevel)
	;


      // Component Wise Residual Convergence Criterion
      class_< ComponentWiseResidualConvergenceCriterionType,
	      bases< ConvergenceCriteriaType >, boost::noncopyable >
	(
	 "ComponentWiseResidualConvergenceCriterion", 
	 init<double, double >())
	.def(init<double, double >())
	.def("SetEchoLevel", &ComponentWiseResidualConvergenceCriterionType::SetEchoLevel)
	;

      
      // Eigensolver Strategy
      class_< EigensolverStrategyType,
	      bases< BaseSolvingStrategyType >, boost::noncopyable >
	(
	 "EigensolverStrategy", init<ModelPart&, BaseSchemeType::Pointer, BuilderAndSolverPointer, bool>() )
	;

    }

  }  // namespace Python.

} // Namespace Kratos

