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
#include "custom_solvers/solution_strategies/newton_raphson_strategy.hpp"
#include "custom_solvers/solution_strategies/line_search_strategy.hpp"
#include "custom_solvers/solution_strategies/explicit_strategy.hpp" 
#include "custom_solvers/solution_strategies/eigensolver_strategy.hpp"

// to update
#include "custom_solvers/solution_strategies/explicit_hamilton_strategy.hpp"

//builders and solvers
#include "custom_solvers/builders_and_solvers/block_builder_and_solver.hpp"
#include "custom_solvers/builders_and_solvers/reduction_builder_and_solver.hpp"
#include "custom_solvers/builders_and_solvers/explicit_builder_and_solver.hpp"

// to update
#include "custom_solvers/builders_and_solvers/explicit_hamilton_builder_and_solver.hpp"

//convergence criteria
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "custom_solvers/convergence_criteria/displacement_convergence_criterion.hpp"

//solution_schemes
#include "custom_solvers/solution_schemes/displacement_bossak_scheme.hpp"
#include "custom_solvers/solution_schemes/eigensolver_dynamic_scheme.hpp" 
 
#include "custom_solvers/solution_schemes/explicit_central_differences_scheme.hpp" 
#include "custom_solvers/solution_schemes/explicit_hamilton_scheme.hpp"

#include "custom_solvers/solution_schemes/displacement_rotation_static_scheme.hpp"
#include "custom_solvers/solution_schemes/displacement_rotation_emc_scheme.hpp"
#include "custom_solvers/solution_schemes/displacement_simo_scheme.hpp"

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
      typedef UblasSpace<double, CompressedMatrix, Vector>                                 SparseSpaceType;
      typedef UblasSpace<double, Matrix, Vector>                                            LocalSpaceType;
      typedef LinearSolver<SparseSpaceType, LocalSpaceType >                              LinearSolverType;
      typedef SolutionStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >   SolutionStrategyType;
      typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType >   BuilderAndSolverType;
      typedef SolutionScheme< SparseSpaceType, LocalSpaceType >                         SolutionSchemeType;
      typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType >               ConvergenceCriteriaType;

      //custom solution strategy types
      typedef LinearStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >                     LinearStrategyType;
      typedef NewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >       NewtonRaphsonStrategyType;
      typedef LineSearchStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >             LineSearchStrategyType;
      typedef ExplicitStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >                 ExplicitStrategyType;
      typedef ExplicitHamiltonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ExplicitHamiltonStrategyType;
      typedef EigensolverStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >           EigensolverStrategyType;

      //custom builder_and_solver types
      typedef ReductionBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType >               ReductionBuilderAndSolverType;
      typedef BlockBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType >                       BlockBuilderAndSolverType;
      typedef ExplicitHamiltonBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType >         ExplicitBuilderAndSolverType;
      typedef ExplicitHamiltonBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > ExplicitHamiltonBuilderAndSolverType;

      //custom solution scheme types
      typedef ExplicitCentralDifferencesScheme< SparseSpaceType, LocalSpaceType >    ExplicitCentralDifferencesSchemeType;
      typedef ExplicitHamiltonScheme< SparseSpaceType, LocalSpaceType >                        ExplicitHamiltonSchemeType;
      typedef EigensolverDynamicScheme< SparseSpaceType, LocalSpaceType >                    EigensolverDynamicSchemeType;

      typedef DisplacementStaticScheme< SparseSpaceType, LocalSpaceType >                    DisplacementStaticSchemeType;     
      typedef DisplacementRotationStaticScheme< SparseSpaceType, LocalSpaceType >    DisplacementRotationStaticSchemeType;
      
      typedef DisplacementNewmarkScheme< SparseSpaceType, LocalSpaceType >                  DisplacementNewmarkSchemeType;     
      typedef DisplacementRotationNewmarkScheme< SparseSpaceType, LocalSpaceType >  DisplacementRotationNewmarkSchemeType;
      
      typedef DisplacementBossakScheme< SparseSpaceType, LocalSpaceType >                    DisplacementBossakSchemeType;
      typedef DisplacementRotationBossakScheme< SparseSpaceType, LocalSpaceType >    DisplacementRotationBossakSchemeType;

      typedef DisplacementSimoScheme< SparseSpaceType, LocalSpaceType >                        DisplacementSimoSchemeType;

      typedef DisplacementRotationSimoScheme< SparseSpaceType, LocalSpaceType >        DisplacementRotationSimoSchemeType;
      typedef DisplacementRotationEmcScheme< SparseSpaceType, LocalSpaceType >          DisplacementRotationEmcSchemeType;
      
      //custom convergence criterion types
      typedef DisplacementConvergenceCriterion< SparseSpaceType,  LocalSpaceType >   DisplacementConvergenceCriterionType;


      //custom integration methods by components
      typedef VariableComponent< VectorComponentAdaptor< array_1d<double, 3 > > >  ComponentType;
      typedef TimeIntegrationMethod<ComponentType, double>                 IntegrationMethodType;
      typedef StaticMethod<ComponentType, double>                               StaticMethodType;
      typedef NewmarkMethod<ComponentType, double>                             NewmarkMethodType;
      typedef BossakMethod<ComponentType, double>                               BossakMethodType;
      typedef SimoMethod<ComponentType, double>                                   SimoMethodType;

      typedef StaticStepMethod<ComponentType, double>                       StaticStepMethodType;	    
      typedef NewmarkStepMethod<ComponentType, double>                     NewmarkStepMethodType;
      typedef BossakStepMethod<ComponentType, double>                       BossakStepMethodType;
      typedef SimoStepMethod<ComponentType, double>                           SimoStepMethodType;
      typedef EmcStepMethod<ComponentType, double>                             EmcStepMethodType;    

      typedef StaticStepRotationMethod<ComponentType, double>       StaticStepRotationMethodType;
      typedef NewmarkStepRotationMethod<ComponentType, double>     NewmarkStepRotationMethodType;
      typedef BossakStepRotationMethod<ComponentType, double>       BossakStepRotationMethodType;
      typedef SimoStepRotationMethod<ComponentType, double>           SimoStepRotationMethodType;
      typedef EmcStepRotationMethod<ComponentType, double>             EmcStepRotationMethodType;
      
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
         .def("GetOptions", &SolutionStrategyType::GetOptions, return_internal_reference<>())
         .def("SetEchoLevel", &SolutionStrategyType::SetEchoLevel)
         .def("GetEchoLevel", &SolutionStrategyType::GetEchoLevel)
         ;

      // Solid Mechanics Linear Strategy
      class_< LinearStrategyType, 
	      bases< SolutionStrategyType >, boost::noncopyable >
          ("LinearStrategy",init < ModelPart&, SolutionSchemeType::Pointer, BuilderAndSolverType::Pointer, Flags& >() )
          .def(init < ModelPart&, SolutionSchemeType::Pointer, LinearSolverType::Pointer, Flags& >())
	;
      
      // Solid Mechanics Newton Raphson Strategy
      class_< NewtonRaphsonStrategyType, 
	      bases< LinearStrategyType >, boost::noncopyable >
          ("NewtonRaphsonStrategy",init < ModelPart&, SolutionSchemeType::Pointer, BuilderAndSolverType::Pointer, ConvergenceCriteriaType::Pointer, Flags&, unsigned int >() )
          .def(init < ModelPart&, SolutionSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, Flags&, unsigned int >())
          .def("SetMaxIterationNumber", &NewtonRaphsonStrategyType::SetMaxIterationNumber)
          .def("GetMaxIterationNumber", &NewtonRaphsonStrategyType::GetMaxIterationNumber)  
	;

      // Solid Mechanics Newton Raphson Line Search Strategy
      class_< LineSearchStrategyType, 
	      bases< NewtonRaphsonStrategyType >, boost::noncopyable >
          ("LineSearchStrategy",init < ModelPart&, SolutionSchemeType::Pointer, BuilderAndSolverType::Pointer, ConvergenceCriteriaType::Pointer, Flags&, unsigned int >() )
          .def(init < ModelPart&, SolutionSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, Flags&, unsigned int >())
	;
      
      // Solid Mechanics Explicit Strategy
      class_< ExplicitStrategyType, 
	      bases< SolutionStrategyType >, boost::noncopyable >
          ("ExplicitStrategy",init < ModelPart&, SolutionSchemeType::Pointer, Flags& >() )
	;
         

      // Eigensolver Strategy
      class_< EigensolverStrategyType,
	      bases< SolutionStrategyType >, boost::noncopyable >
	(
            "EigensolverStrategy", init<ModelPart&, SolutionSchemeType::Pointer, BuilderAndSolverType::Pointer, Flags&, bool>() )
	;

      
      // // Explicit Hamilton Estrategy for Explicit Beam solution
      // class_< ExplicitHamiltonStrategyType, 
      //         bases< BaseSolvingStrategyType >, boost::noncopyable >
      //   (
      //    "ExplicitHamiltonStrategy",
      //    init < ModelPart&, SolutionSchemeType::Pointer,  LinearSolverType::Pointer, bool, bool, bool>())
	
      //   .def(init < ModelPart&, SolutionSchemeType::Pointer, LinearSolverType::Pointer,  bool, bool, bool >())
      //   ;

      
      //********************************************************************
      //*******************BUILDER AND SOLVER CLASSES***********************
      //********************************************************************

      // Solid Mechanics Base Builder and Solver
      class_< BuilderAndSolverType, bases<Flags>, boost::noncopyable >
          ("BuilderAndSolver", init<LinearSolverType::Pointer > ())
          .def(init<>())
          .def("BuildLHS", &BuilderAndSolverType::BuildLHS)
          .def("BuildRHS", &BuilderAndSolverType::BuildRHS)
          .def("Build", &BuilderAndSolverType::Build)
          .def("SystemSolve", &BuilderAndSolverType::SystemSolve)
          .def("BuildAndSolve", &BuilderAndSolverType::BuildAndSolve)
          .def("BuildRHSAndSolve", &BuilderAndSolverType::BuildRHSAndSolve)
          .def("SetUpDofSet", &BuilderAndSolverType::SetUpDofSet)
          .def("GetDofSet", &BuilderAndSolverType::GetDofSet, return_internal_reference<>())
          .def("SetUpSystem", &BuilderAndSolverType::SetUpSystem)
          .def("SetUpSystemMatrices", &BuilderAndSolverType::SetUpSystemMatrices)
          .def("InitializeSolutionStep", &BuilderAndSolverType::InitializeSolutionStep)
          .def("FinalizeSolutionStep", &BuilderAndSolverType::FinalizeSolutionStep)
          .def("CalculateReactions", &BuilderAndSolverType::CalculateReactions)
          .def("GetEquationSystemSize", &BuilderAndSolverType::GetEquationSystemSize)
          .def("Clear", &BuilderAndSolverType::Clear)
          .def("Check", &BuilderAndSolverType::Check)
          .def("SetEchoLevel", &BuilderAndSolverType::SetEchoLevel)
          .def("GetEchoLevel", &BuilderAndSolverType::GetEchoLevel)
          ;

      class_< ReductionBuilderAndSolverType, bases<BuilderAndSolverType>, boost::noncopyable >
        (
         "ReductionBuilderAndSolver", init< LinearSolverType::Pointer > ()
         );
            
      class_< BlockBuilderAndSolverType, bases<BuilderAndSolverType>, boost::noncopyable >
        (
         "BlockBuilderAndSolver", init< LinearSolverType::Pointer > ()
         );

      class_< ExplicitBuilderAndSolverType, bases<BuilderAndSolverType>, boost::noncopyable > 
	(
	 "ExplicitBuilderAndSolver", init<> ()
	 );
      
 
      // class_< ExplicitHamiltonBuilderAndSolverType, bases<BuilderAndSolverType>, boost::noncopyable > 
      //   (
      //    "ExplicitHamiltonBuilderAndSolver", init< LinearSolverType::Pointer > ()
      //    );

      
      //********************************************************************
      //*************************SHCHEME CLASSES****************************
      //********************************************************************

      // Solid Mechanics Base Scheme      
      class_< SolutionSchemeType, bases<Flags>, boost::noncopyable >
          ("SolutionScheme", init< >())
          .def("Initialize", &SolutionSchemeType::Initialize)
          .def("InitializeElements", &SolutionSchemeType::InitializeElements)
          .def("InitializeConditions", &SolutionSchemeType::InitializeConditions)
          .def("InitializeSolutionStep", &SolutionSchemeType::InitializeSolutionStep)
          .def("FinalizeSolutionStep", &SolutionSchemeType::FinalizeSolutionStep)
          .def("InitializeNonLinearIteration", &SolutionSchemeType::InitializeNonLinearIteration)
          .def("FinalizeNonLinearIteration", &SolutionSchemeType::FinalizeNonLinearIteration)
          .def("Predict", &SolutionSchemeType::Predict)
          .def("Update", &SolutionSchemeType::Update)
          .def("MoveMesh", &SolutionSchemeType::MoveMesh)
          .def("Clear",&SolutionSchemeType::Clear)
          .def("Check", &SolutionSchemeType::Check)
          ;
      
      // Explicit scheme: Central differences 
      class_< ExplicitCentralDifferencesSchemeType,
	      bases< SolutionSchemeType >,  boost::noncopyable >
	(
	 "ExplicitCentralDifferencesScheme", init< const double, const double, const double, const bool >() )
	;

      // // Explicit Hamilton Scheme Type
      // class_< ExplicitHamiltonSchemeType,
      //         bases< SolutionSchemeType >,  boost::noncopyable >
      //   (
      //    "ExplicitHamiltonScheme", init< double, double, double, bool >() )
	
      //   ;


      // Displacement Static Scheme Type
      class_< DisplacementStaticSchemeType,
      	      bases< SolutionSchemeType >,  boost::noncopyable >
      	(
      	 "DisplacementStaticScheme", init<>() )
      	;
	    
      // Displacement Rotation Static Scheme Type
      class_< DisplacementRotationStaticSchemeType,
      	      bases< SolutionSchemeType >,  boost::noncopyable >
      	(
      	 "DisplacementRotationStaticScheme", init<>() )
      	;

      // Displacement Newmark Scheme Type
      class_< DisplacementNewmarkSchemeType,
      	      bases< SolutionSchemeType >,  boost::noncopyable >
      	(
      	 "DisplacementNewmarkScheme", init<>() )
      	;
	    
      // Displacement Rotation Newmark Scheme Type
      class_< DisplacementRotationNewmarkSchemeType,
      	      bases< SolutionSchemeType >,  boost::noncopyable >
      	(
      	 "DisplacementRotationNewmarkScheme", init<>() )
      	;
      
      // Displacement Bossak Scheme Type
      class_< DisplacementBossakSchemeType,
      	      bases< SolutionSchemeType >,  boost::noncopyable >
      	(
      	 "DisplacementBossakScheme", init<>() )
      	;
      
      // Displacement Rotation Bossak Scheme Type
      class_< DisplacementRotationBossakSchemeType,
      	      bases< SolutionSchemeType >,  boost::noncopyable >
      	(
      	 "DisplacementRotationBossakScheme", init<>() )
      	;

      // Displacement Simo Scheme Type
      class_< DisplacementSimoSchemeType,
      	      bases< SolutionSchemeType >,  boost::noncopyable >
      	(
      	 "DisplacementSimoScheme", init<>() )
       	;

      // Displacement Rotation Simo Scheme Type
      class_< DisplacementRotationSimoSchemeType,
      	      bases< SolutionSchemeType >,  boost::noncopyable >
      	(
      	 "DisplacementRotationSimoScheme", init<>() )
      	;
      
      // Displacement Rotation Emc Scheme Type
      class_< DisplacementRotationEmcSchemeType,
      	      bases< SolutionSchemeType >,  boost::noncopyable >
      	(
      	 "DisplacementRotationEmcScheme", init<>() )
      	;


      // Eigensolver Scheme Type
      class_< EigensolverDynamicSchemeType,
	      bases< SolutionSchemeType >, boost::noncopyable >
	(
	 "EigensolverDynamicScheme", init<>() )
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

