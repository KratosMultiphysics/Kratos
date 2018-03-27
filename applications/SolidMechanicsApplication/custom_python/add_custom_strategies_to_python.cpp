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
#include "custom_strategies/strategies/linear_strategy.hpp"
#include "custom_strategies/strategies/newton_raphson_strategy.hpp"
#include "custom_strategies/strategies/component_wise_newton_raphson_strategy.hpp"
#include "custom_strategies/strategies/line_search_strategy.hpp"
#include "custom_strategies/strategies/explicit_strategy.hpp" 
#include "custom_strategies/strategies/eigensolver_strategy.hpp"

// to update
#include "custom_strategies/strategies/explicit_hamilton_strategy.hpp"

//builders and solvers
#include "custom_strategies/builders_and_solvers/component_wise_builder_and_solver.hpp"
#include "custom_strategies/builders_and_solvers/explicit_hamilton_builder_and_solver.hpp"

//convergence criteria
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "custom_strategies/convergence_criteria/displacement_convergence_criterion.hpp"
#include "custom_strategies/convergence_criteria/component_wise_residual_convergence_criterion.hpp"

//schemes
#include "custom_strategies/schemes/residual_based_displacement_bossak_scheme.hpp"
#include "custom_strategies/schemes/component_wise_bossak_scheme.hpp"
#include "custom_strategies/schemes/eigensolver_dynamic_scheme.hpp" 
 
#include "custom_strategies/schemes/explicit_central_differences_scheme.hpp" 
#include "custom_strategies/schemes/residual_based_rotation_newmark_scheme.hpp"
#include "custom_strategies/schemes/residual_based_rotation_simo_scheme.hpp"
#include "custom_strategies/schemes/residual_based_rotation_emc_scheme.hpp"
#include "custom_strategies/schemes/explicit_hamilton_scheme.hpp"

#include "custom_strategies/schemes/residual_based_displacement_rotation_static_scheme.hpp"
#include "custom_strategies/schemes/residual_based_displacement_rotation_emc_scheme.hpp"
#include "custom_strategies/schemes/residual_based_displacement_simo_scheme.hpp"

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
      //typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
      typedef SolutionStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > SolutionStrategyType;
      typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;
      typedef Scheme< SparseSpaceType, LocalSpaceType > SchemeType;
      typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaType;

      //custom strategy types
      typedef LinearStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > LinearStrategyType;
      typedef NewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > NewtonRaphsonStrategyType;
      typedef ComponentWiseNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ComponentWiseNewtonRaphsonStrategyType;
      typedef LineSearchStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > LineSearchStrategyType;
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

      typedef ResidualBasedDisplacementStaticScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedDisplacementStaticSchemeType;     
      typedef ResidualBasedDisplacementRotationStaticScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedDisplacementRotationStaticSchemeType;
      
      typedef ResidualBasedDisplacementNewmarkScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedDisplacementNewmarkSchemeType;     
      typedef ResidualBasedDisplacementRotationNewmarkScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedDisplacementRotationNewmarkSchemeType;
      
      typedef ResidualBasedDisplacementBossakScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedDisplacementBossakSchemeType;
      typedef ResidualBasedDisplacementRotationBossakScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedDisplacementRotationBossakSchemeType;

      typedef ResidualBasedDisplacementSimoScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedDisplacementSimoSchemeType;

      typedef ResidualBasedDisplacementRotationSimoScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedDisplacementRotationSimoSchemeType;
      typedef ResidualBasedDisplacementRotationEmcScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedDisplacementRotationEmcSchemeType;
      
      //custom convergence criterion types
      typedef DisplacementConvergenceCriterion< SparseSpaceType,  LocalSpaceType > DisplacementConvergenceCriterionType;
      typedef ComponentWiseResidualConvergenceCriterion< SparseSpaceType,  LocalSpaceType > ComponentWiseResidualConvergenceCriterionType;


      //custom integration methods by components
      typedef VariableComponent< VectorComponentAdaptor< array_1d<double, 3 > > >  ComponentType;
      typedef TimeIntegrationMethod<ComponentType, double>  IntegrationMethodType;
      typedef StaticMethod<ComponentType, double>  StaticMethodType;
      typedef NewmarkMethod<ComponentType, double>  NewmarkMethodType;
      typedef BossakMethod<ComponentType, double>  BossakMethodType;
      typedef SimoMethod<ComponentType, double>  SimoMethodType;

      typedef StaticStepMethod<ComponentType, double>  StaticStepMethodType;	    
      typedef NewmarkStepMethod<ComponentType, double>  NewmarkStepMethodType;
      typedef BossakStepMethod<ComponentType, double>  BossakStepMethodType;
      typedef SimoStepMethod<ComponentType, double>  SimoStepMethodType;
      typedef EmcStepMethod<ComponentType, double>  EmcStepMethodType;    

      typedef StaticStepRotationMethod<ComponentType, double>  StaticStepRotationMethodType;
      typedef NewmarkStepRotationMethod<ComponentType, double>  NewmarkStepRotationMethodType;
      typedef BossakStepRotationMethod<ComponentType, double>  BossakStepRotationMethodType;
      typedef SimoStepRotationMethod<ComponentType, double>  SimoStepRotationMethodType;
      typedef EmcStepRotationMethod<ComponentType, double>  EmcStepRotationMethodType;
      
      //********************************************************************
      //*************************STRATEGY CLASSES***************************
      //********************************************************************

      // Solver Local Flags
      class_<SolverLocalFlags, boost::noncopyable> ("SolverLocalFlags", init<>() )
          .def_readonly("INITIALIZED", &SolverLocalFlags::INITIALIZED)
          .def_readonly("MOVE_MESH", &SolverLocalFlags::MOVE_MESH)
          .def_readonly("REFORM_DOFS", &SolverLocalFlags::REFORM_DOFS)
          .def_readonly("CONSTANT_SYSTEM_LHS", &SolverLocalFlags::CONSTANT_SYSTEM_LHS)
          .def_readonly("COMPUTE_REACTIONS", &SolverLocalFlags::COMPUTE_REACTIONS)
          .def_readonly("IMPLEX", &SolverLocalFlags::IMPLEX)
          ;
      
      // Solid Mechanics Base Strategy
      class_< SolutionStrategyType, bases<Flags>, boost::noncopyable >("SolutionStrategy", init< ModelPart& >() )
         .def(init < ModelPart&, Flags& >())
         .def("InitializeSolutionStep", &SolutionStrategyType::InitializeSolutionStep)
         .def("FinalizeSolutionStep", &SolutionStrategyType::FinalizeSolutionStep)
         .def("SolveSolutionStep", &SolutionStrategyType::SolveSolutionStep )
         .def("Solve", &SolutionStrategyType::Solve)
         .def("Check", &SolutionStrategyType::Check)
         .def("Clear", &SolutionStrategyType::Clear)
         .def("SetOptions", &SolutionStrategyType::SetOptions)
         .def("SetEchoLevel", &SolutionStrategyType::SetEchoLevel)
         .def("GetEchoLevel", &SolutionStrategyType::GetEchoLevel)
         ;

      // Solid Mechanics Linear Strategy
      class_< LinearStrategyType, 
	      bases< SolutionStrategyType >, boost::noncopyable >
          ("LinearStrategy",init < ModelPart&, SchemeType::Pointer, BuilderAndSolverType::Pointer, Flags& >() )
          .def(init < ModelPart&, SchemeType::Pointer, LinearSolverType::Pointer, Flags& >())
	;
      
      // Solid Mechanics Newton Raphson Strategy
      class_< NewtonRaphsonStrategyType, 
	      bases< SolutionStrategyType >, boost::noncopyable >
          ("NewtonRaphsonStrategy",init < ModelPart&, SchemeType::Pointer, BuilderAndSolverType::Pointer, ConvergenceCriteriaType::Pointer, Flags&, unsigned int >() )
          .def(init < ModelPart&, SchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, Flags&, unsigned int >())
          .def("SetMaxIterationNumber", &NewtonRaphsonStrategyType::SetMaxIterationNumber)
          .def("GetMaxIterationNumber", &NewtonRaphsonStrategyType::GetMaxIterationNumber)  
	;

      // Solid Mechanics Newton Raphson Line Search Strategy
      class_< LineSearchStrategyType, 
	      bases< NewtonRaphsonStrategyType >, boost::noncopyable >
          ("LineSearchStrategy",init < ModelPart&, SchemeType::Pointer, BuilderAndSolverType::Pointer, ConvergenceCriteriaType::Pointer, Flags&, unsigned int >() )
          .def(init < ModelPart&, SchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, Flags&, unsigned int >())
	;
      
      // Solid Mechanics Explicit Strategy
      class_< ExplicitStrategyType, 
	      bases< SolutionStrategyType >, boost::noncopyable >
          ("ExplicitStrategy",init < ModelPart&, SchemeType::Pointer, Flags& >() )
	;
        
      // Component Wise Newton-Raphson Strategy
      class_< ComponentWiseNewtonRaphsonStrategyType, 
	      bases< NewtonRaphsonStrategyType >, boost::noncopyable >
	(
	 "ComponentWiseNewtonRaphsonStrategy",
	 init <ModelPart&, SchemeType::Pointer, BuilderAndSolverType::Pointer, ConvergenceCriteriaType::Pointer, Flags&, unsigned int >() )
          .def(init < ModelPart&, SchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, Flags&, unsigned int >())
	;
 

      // Eigensolver Strategy
      class_< EigensolverStrategyType,
	      bases< SolutionStrategyType >, boost::noncopyable >
	(
            "EigensolverStrategy", init<ModelPart&, SchemeType::Pointer, BuilderAndSolverType::Pointer, Flags&, bool>() )
	;

      
      // // Explicit Hamilton Estrategy for Explicit Beam solution
      // class_< ExplicitHamiltonStrategyType, 
      //         bases< BaseSolvingStrategyType >, boost::noncopyable >
      //   (
      //    "ExplicitHamiltonStrategy",
      //    init < ModelPart&, SchemeType::Pointer,  LinearSolverType::Pointer, bool, bool, bool>())
	
      //   .def(init < ModelPart&, SchemeType::Pointer, LinearSolverType::Pointer,  bool, bool, bool >())
      //   ;

      
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
	      bases< SchemeType >,  boost::noncopyable >
	(
	 "ComponentWiseBossakScheme", init< double >() )

	.def("Initialize", &ComponentWiseBossakScheme<SparseSpaceType, LocalSpaceType>::Initialize)
	;


      // Explicit scheme: Central differences 
      class_< ExplicitCentralDifferencesSchemeType,
	      bases< SchemeType >,  boost::noncopyable >
	(
	 "ExplicitCentralDifferencesScheme", init< const double, const double, const double, const bool >() )

	.def("Initialize", &ExplicitCentralDifferencesScheme<SparseSpaceType, LocalSpaceType>::Initialize)
	;

      // Residual Based Rotational Newmark Scheme Type
      class_< ResidualBasedRotationNewmarkSchemeType,
	      bases< SchemeType >,  boost::noncopyable >
	(
	 "ResidualBasedRotationNewmarkScheme", init< double, double >() )
	
	.def("Initialize", &ResidualBasedRotationNewmarkScheme<SparseSpaceType, LocalSpaceType>::Initialize)
	;

      // Residual Based Rotational Simo Scheme Type
      class_< ResidualBasedRotationSimoSchemeType,
	      bases< SchemeType >,  boost::noncopyable >
	(
	 "ResidualBasedRotationSimoScheme", init< double, double >() )
	
	.def("Initialize", &ResidualBasedRotationSimoScheme<SparseSpaceType, LocalSpaceType>::Initialize)
	;

      // Residual Based Rotational EMC Scheme Type
      class_< ResidualBasedRotationEMCSchemeType,
	      bases< SchemeType >,  boost::noncopyable >
	(
	 "ResidualBasedRotationEMCScheme", init< double >() )
	
	.def("Initialize", &ResidualBasedRotationEMCScheme<SparseSpaceType, LocalSpaceType>::Initialize)
	;

      // Explicit Hamilton Scheme Type
      class_< ExplicitHamiltonSchemeType,
	      bases< SchemeType >,  boost::noncopyable >
	(
	 "ExplicitHamiltonScheme", init< double, double, double, bool >() )
	
	.def("Initialize", &ExplicitHamiltonScheme<SparseSpaceType, LocalSpaceType>::Initialize)
	;


      // Residual Based Displacement Static Scheme Type
      class_< ResidualBasedDisplacementStaticSchemeType,
      	      bases< SchemeType >,  boost::noncopyable >
      	(
      	 "ResidualBasedDisplacementStaticScheme", init<>() )
      	.def("Initialize", &ResidualBasedDisplacementStaticScheme<SparseSpaceType, LocalSpaceType>::Initialize)
      	;
	    
      // Residual Based Displacement Rotation Static Scheme Type
      class_< ResidualBasedDisplacementRotationStaticSchemeType,
      	      bases< SchemeType >,  boost::noncopyable >
      	(
      	 "ResidualBasedDisplacementRotationStaticScheme", init<>() )
      	.def("Initialize", &ResidualBasedDisplacementRotationStaticScheme<SparseSpaceType, LocalSpaceType>::Initialize)
      	;

      // Residual Based Displacement Newmark Scheme Type
      class_< ResidualBasedDisplacementNewmarkSchemeType,
      	      bases< SchemeType >,  boost::noncopyable >
      	(
      	 "ResidualBasedDisplacementNewmarkScheme", init<>() )
      	.def("Initialize", &ResidualBasedDisplacementNewmarkScheme<SparseSpaceType, LocalSpaceType>::Initialize)
      	;
	    
      // Residual Based Displacement Rotation Newmark Scheme Type
      class_< ResidualBasedDisplacementRotationNewmarkSchemeType,
      	      bases< SchemeType >,  boost::noncopyable >
      	(
      	 "ResidualBasedDisplacementRotationNewmarkScheme", init<>() )
      	.def("Initialize", &ResidualBasedDisplacementRotationNewmarkScheme<SparseSpaceType, LocalSpaceType>::Initialize)
      	;
      
      // Residual Based Displacement Bossak Scheme Type
      class_< ResidualBasedDisplacementBossakSchemeType,
      	      bases< SchemeType >,  boost::noncopyable >
      	(
      	 "ResidualBasedDisplacementBossakScheme", init<>() )
      	.def("Initialize", &ResidualBasedDisplacementBossakScheme<SparseSpaceType, LocalSpaceType>::Initialize)
      	;
      
      // Residual Based Displacement Rotation Bossak Scheme Type
      class_< ResidualBasedDisplacementRotationBossakSchemeType,
      	      bases< SchemeType >,  boost::noncopyable >
      	(
      	 "ResidualBasedDisplacementRotationBossakScheme", init<>() )
      	.def("Initialize", &ResidualBasedDisplacementRotationBossakScheme<SparseSpaceType, LocalSpaceType>::Initialize)
      	;

      // Residual Based Displacement Simo Scheme Type
      class_< ResidualBasedDisplacementSimoSchemeType,
      	      bases< SchemeType >,  boost::noncopyable >
      	(
      	 "ResidualBasedDisplacementSimoScheme", init<>() )
      	.def("Initialize", &ResidualBasedDisplacementSimoScheme<SparseSpaceType, LocalSpaceType>::Initialize)
      	;

      // Residual Based Displacement Rotation Simo Scheme Type
      class_< ResidualBasedDisplacementRotationSimoSchemeType,
      	      bases< SchemeType >,  boost::noncopyable >
      	(
      	 "ResidualBasedDisplacementRotationSimoScheme", init<>() )
      	.def("Initialize", &ResidualBasedDisplacementRotationSimoScheme<SparseSpaceType, LocalSpaceType>::Initialize)
      	;
      
      // Residual Based Displacement Rotation Emc Scheme Type
      class_< ResidualBasedDisplacementRotationEmcSchemeType,
      	      bases< SchemeType >,  boost::noncopyable >
      	(
      	 "ResidualBasedDisplacementRotationEmcScheme", init<>() )
      	.def("Initialize", &ResidualBasedDisplacementRotationEmcScheme<SparseSpaceType, LocalSpaceType>::Initialize)
      	;


      // Eigensolver Scheme Type
      class_< EigensolverDynamicSchemeType, EigensolverDynamicSchemeType::Pointer,
	      bases< SchemeType >, boost::noncopyable >
	(
	 "EigensolverDynamicScheme",
	 init<>())
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
	.def("SetEchoLevel", &DisplacementConvergenceCriterionType::SetEchoLevel)
	;


      // Component Wise Residual Convergence Criterion
      class_< ComponentWiseResidualConvergenceCriterionType,
	      bases< ConvergenceCriteriaType >, boost::noncopyable >
	(
	 "ComponentWiseResidualConvergenceCriterion", 
	 init<double, double >())
	.def("SetEchoLevel", &ComponentWiseResidualConvergenceCriterionType::SetEchoLevel)
	;

      


      //********************************************************************
      //*******************TIME INTEGRATION METHODS*************************
      //********************************************************************

      //Time Integration Method   
      class_< IntegrationMethodType, IntegrationMethodType::Pointer, boost::noncopyable >
      	(
      	 "TimeIntegrationMethod", init<>())
	.def("Clone", &IntegrationMethodType::Clone)
	.def("SetVariable", &IntegrationMethodType::SetVariable)
	.def("SetFirstDerivative", &IntegrationMethodType::SetFirstDerivative)
	.def("SetSecondDerivative", &IntegrationMethodType::SetSecondDerivative)
      	.def("SetVariables", &IntegrationMethodType::SetVariables)
      	.def("SetInputVariable", &IntegrationMethodType::SetInputVariable)
	.def("HasStepVariable", &IntegrationMethodType::HasStepVariable)
	.def("SetStepVariable", &IntegrationMethodType::SetStepVariable)
        .def("CalculateParameters", &IntegrationMethodType::CalculateParameters)
      	.def("SetParameters", &IntegrationMethodType::SetParameters)
      	.def("Predict", &IntegrationMethodType::Predict)
	.def(self_ns::str(self))
	DECLARE_HAS_THIS_TYPE_PROCESS_INFO_PYTHON_AS_POINTER(IntegrationMethodType)
	DECLARE_ADD_THIS_TYPE_TO_PROCESS_INFO_PYTHON_AS_POINTER(IntegrationMethodType)
	DECLARE_GET_THIS_TYPE_FROM_PROCESS_INFO_PYTHON_AS_POINTER(IntegrationMethodType)
      	;

      //to define it as a variable 
      class_<Variable<IntegrationMethodType::Pointer> , bases<VariableData>, boost::noncopyable >("TimeIntegrationMethodVariable", no_init)
	;

      class_< StaticMethodType, StaticMethodType::Pointer,
      	      bases< IntegrationMethodType >, boost::noncopyable >
      	(
      	 "StaticMethod", init<>())
       	;
      
      class_< NewmarkMethodType, NewmarkMethodType::Pointer,
      	      bases< IntegrationMethodType >, boost::noncopyable >
      	(
      	 "NewmarkMethod", init<>())
       	;

      class_< BossakMethodType, BossakMethodType::Pointer,
      	      bases< IntegrationMethodType >, boost::noncopyable >
      	(
      	 "BossakMethod", init<>())
      	;
      
      class_< SimoMethodType, SimoMethodType::Pointer,
      	      bases< IntegrationMethodType >, boost::noncopyable >
      	(
      	 "SimoMethod", init<>())
      	;
      
      class_< StaticStepMethodType, StaticStepMethodType::Pointer,
      	      bases< IntegrationMethodType >, boost::noncopyable >
      	(
      	 "StaticStepMethod", init<>())
     	;

      class_< NewmarkStepMethodType, NewmarkStepMethodType::Pointer,
      	      bases< IntegrationMethodType >, boost::noncopyable >
      	(
      	 "NewmarkStepMethod", init<>())
      	;

      class_< BossakStepMethodType, BossakStepMethodType::Pointer,
      	      bases< IntegrationMethodType >, boost::noncopyable >
      	(
      	 "BossakStepMethod", init<>())
      	;

      class_< SimoStepMethodType, SimoStepMethodType::Pointer,
      	      bases< IntegrationMethodType >, boost::noncopyable >
      	(
      	 "SimoStepMethod", init<>())
      	;
      
       class_< EmcStepMethodType, EmcStepMethodType::Pointer,
      	      bases< IntegrationMethodType >, boost::noncopyable >
      	(
      	 "EmcStepMethod", init<>())
      	;

       class_< StaticStepRotationMethodType, StaticStepRotationMethodType::Pointer,
      	      bases< IntegrationMethodType >, boost::noncopyable >
      	(
      	 "StaticStepRotationMethod", init<>())
      	;

       class_< NewmarkStepRotationMethodType, NewmarkStepRotationMethodType::Pointer,
      	      bases< IntegrationMethodType >, boost::noncopyable >
      	(
      	 "NewmarkStepRotationMethod", init<>())
      	;

       class_< BossakStepRotationMethodType, BossakStepRotationMethodType::Pointer,
      	      bases< IntegrationMethodType >, boost::noncopyable >
      	(
      	 "BossakStepRotationMethod", init<>())
      	;

       class_< SimoStepRotationMethodType, SimoStepRotationMethodType::Pointer,
      	      bases< IntegrationMethodType >, boost::noncopyable >
      	(
      	 "SimoStepRotationMethod", init<>())
       	;

       class_< EmcStepRotationMethodType, EmcStepRotationMethodType::Pointer,
      	      bases< IntegrationMethodType >, boost::noncopyable >
      	(
      	 "EmcStepRotationMethod", init<>())
      	;
    }

  }  // namespace Python.

} // Namespace Kratos

