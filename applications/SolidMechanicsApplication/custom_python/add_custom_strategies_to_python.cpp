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
#include "custom_solvers/solution_strategies/segregated_strategy.hpp"

// to update
#include "custom_solvers/solution_strategies/explicit_hamilton_strategy.hpp"

//builders and solvers
#include "custom_solvers/solution_builders_and_solvers/block_builder_and_solver.hpp"
#include "custom_solvers/solution_builders_and_solvers/reduction_builder_and_solver.hpp"
#include "custom_solvers/solution_builders_and_solvers/explicit_builder_and_solver.hpp"

// to update
#include "custom_solvers/solution_builders_and_solvers/explicit_hamilton_builder_and_solver.hpp"

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

#include "custom_solvers/time_integration_methods/time_integration_methods_container.hpp"


namespace Kratos
{

  namespace Python
  {
    using namespace boost::python;

    //base types
    typedef UblasSpace<double, CompressedMatrix, Vector>                                               SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector>                                                          LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType >                                            LinearSolverType;
    typedef SolutionStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >                 SolutionStrategyType;
    typedef SolutionBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > SolutionBuilderAndSolverType;
    typedef SolutionScheme< SparseSpaceType, LocalSpaceType >                                       SolutionSchemeType;
    typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType >                             ConvergenceCriteriaType; 
    typedef SolutionStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >                 SolutionStrategyType;
  
    typedef SolutionStrategyType::Pointer                                                      SolutionStrategyPointer;
    typedef std::vector<SolutionStrategyType::Pointer>                                     SolutionStrategiesContainer;
  
    void Push_Back_Solution_Strategies( SolutionStrategiesContainer& ThisSolutionStrategyContainer,
                                        SolutionStrategyPointer ThisSolutionStrategy )
    { 
      ThisSolutionStrategyContainer.push_back( ThisSolutionStrategy );
    } 

    void  AddCustomStrategiesToPython()
    {

      //custom solution strategy types
      typedef SegregatedStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >             SegregatedStrategyType;
      typedef LinearStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >                     LinearStrategyType;
      typedef NewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >       NewtonRaphsonStrategyType;
      typedef LineSearchStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >             LineSearchStrategyType;
      typedef ExplicitStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >                 ExplicitStrategyType;
      typedef ExplicitHamiltonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ExplicitHamiltonStrategyType;
      typedef EigensolverStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >           EigensolverStrategyType;

      //custom builder_and_solver types
      typedef ReductionBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType >               ReductionBuilderAndSolverType;
      typedef BlockBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType >                       BlockBuilderAndSolverType;
      typedef ExplicitBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType >                 ExplicitBuilderAndSolverType;
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

      //custom integration methods by variables
      typedef array_1d<double, 3>                                                            VectorType;
      typedef Variable<VectorType >                                                        VariableType; 
      typedef TimeIntegrationMethod<VariableType, VectorType>                 TimeIntegrationMethodType;
      typedef typename TimeIntegrationMethodType::Pointer              TimeIntegrationMethodPointerType;

      //custom integration methods by components
      typedef TimeIntegrationMethodsContainer                              TimeIntegrationContainerType;      
      typedef TimeIntegrationContainerType::Pointer                 TimeIntegrationContainerPointerType;
      
      typedef VariableComponent< VectorComponentAdaptor< VectorType > >           VariableComponentType;
      typedef TimeIntegrationMethod<VariableComponentType, double>   TimeIntegrationComponentMethodType;      
      typedef StaticMethod<VariableComponentType, double>                              StaticMethodType;
      typedef NewmarkMethod<VariableComponentType, double>                            NewmarkMethodType;
      typedef BossakMethod<VariableComponentType, double>                              BossakMethodType;
      typedef SimoMethod<VariableComponentType, double>                                  SimoMethodType;

      typedef StaticStepMethod<VariableComponentType, double>                      StaticStepMethodType;	    
      typedef NewmarkStepMethod<VariableComponentType, double>                    NewmarkStepMethodType;
      typedef BossakStepMethod<VariableComponentType, double>                      BossakStepMethodType;
      typedef SimoStepMethod<VariableComponentType, double>                          SimoStepMethodType;
      typedef EmcStepMethod<VariableComponentType, double>                            EmcStepMethodType;    

      typedef StaticStepRotationMethod<VariableComponentType, double>      StaticStepRotationMethodType;
      typedef NewmarkStepRotationMethod<VariableComponentType, double>    NewmarkStepRotationMethodType;
      typedef BossakStepRotationMethod<VariableComponentType, double>      BossakStepRotationMethodType;
      typedef SimoStepRotationMethod<VariableComponentType, double>          SimoStepRotationMethodType;
      typedef EmcStepRotationMethod<VariableComponentType, double>            EmcStepRotationMethodType;

     
      //***************************SOLVER FLAGS******************************

      // Solver Local Flags
      class_<SolverLocalFlags, boost::noncopyable> ("SolverLocalFlags", init<>() )
          .def_readonly("INITIALIZED", &SolverLocalFlags::INITIALIZED)
          .def_readonly("CONVERGED", &SolverLocalFlags::CONVERGED)
          .def_readonly("REFORM_DOFS", &SolverLocalFlags::REFORM_DOFS)
          .def_readonly("COMPUTE_REACTIONS", &SolverLocalFlags::COMPUTE_REACTIONS)
          .def_readonly("CONSTANT_SYSTEM_MATRIX", &SolverLocalFlags::CONSTANT_SYSTEM_MATRIX)
          .def_readonly("RAYLEIGH_DAMPING", &SolverLocalFlags::RAYLEIGH_DAMPING)
          .def_readonly("IMPLEX", &SolverLocalFlags::IMPLEX)
          ;

      
      //*************************STRATEGY CLASSES***************************

      // Solid Mechanics Solution Strategies Container
      class_< SolutionStrategiesContainer >( "SolutionStragetiesContainer", init<>() )
          .def( "PushBack", Push_Back_Solution_Strategies )
          ;
      
      // Solid Mechanics Base Solution Strategy
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

      // Solid Mechanics Segregated Strategy
      class_< SegregatedStrategyType, 
	      bases< SolutionStrategyType >, boost::noncopyable >
          ("SegregatedStrategy",init < ModelPart&, Flags& >() )
          .def(init< ModelPart&, Flags&,  SolutionStrategiesContainer& >())
          .def("AddStrategy", &SegregatedStrategyType::AddStrategy)
          .def("GetStrategy", &SegregatedStrategyType::GetStrategy)
	;
      
      // Solid Mechanics Linear Strategy
      class_< LinearStrategyType, 
	      bases< SolutionStrategyType >, boost::noncopyable >
          ("LinearStrategy",init < ModelPart&, SolutionSchemeType::Pointer, SolutionBuilderAndSolverType::Pointer, Flags& >() )
          .def(init < ModelPart&, SolutionSchemeType::Pointer, LinearSolverType::Pointer, Flags& >())
	;
      
      // Solid Mechanics Newton Raphson Strategy
      class_< NewtonRaphsonStrategyType, 
	      bases< LinearStrategyType >, boost::noncopyable >
          ("NewtonRaphsonStrategy",init < ModelPart&, SolutionSchemeType::Pointer, SolutionBuilderAndSolverType::Pointer, ConvergenceCriteriaType::Pointer, Flags&, unsigned int >() )
          .def(init < ModelPart&, SolutionSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, Flags&, unsigned int >())
          .def("SetMaxIterationNumber", &NewtonRaphsonStrategyType::SetMaxIterationNumber)
          .def("GetMaxIterationNumber", &NewtonRaphsonStrategyType::GetMaxIterationNumber)  
	;

      // Solid Mechanics Newton Raphson Line Search Strategy
      class_< LineSearchStrategyType, 
	      bases< NewtonRaphsonStrategyType >, boost::noncopyable >
          ("LineSearchStrategy",init < ModelPart&, SolutionSchemeType::Pointer, SolutionBuilderAndSolverType::Pointer, ConvergenceCriteriaType::Pointer, Flags&, unsigned int >() )
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
            "EigensolverStrategy", init<ModelPart&, SolutionSchemeType::Pointer, SolutionBuilderAndSolverType::Pointer, Flags&, bool>() )
	;

      
      // // Explicit Hamilton Estrategy for Explicit Beam solution
      // class_< ExplicitHamiltonStrategyType, 
      //         bases< BaseSolvingStrategyType >, boost::noncopyable >
      //   (
      //    "ExplicitHamiltonStrategy",
      //    init < ModelPart&, SolutionSchemeType::Pointer,  LinearSolverType::Pointer, bool, bool, bool>())
	
      //   .def(init < ModelPart&, SolutionSchemeType::Pointer, LinearSolverType::Pointer,  bool, bool, bool >())
      //   ;

      
      //*******************BUILDER AND SOLVER CLASSES***********************

      // Solid Mechanics Base Builder and Solver
      class_< SolutionBuilderAndSolverType, bases<Flags>, boost::noncopyable >
          ("SolutionBuilderAndSolver", init<LinearSolverType::Pointer > ())
          .def(init<>())
          .def("BuildLHS", &SolutionBuilderAndSolverType::BuildLHS)
          .def("BuildRHS", &SolutionBuilderAndSolverType::BuildRHS)
          .def("Build", &SolutionBuilderAndSolverType::Build)
          .def("SystemSolve", &SolutionBuilderAndSolverType::SystemSolve)
          .def("BuildAndSolve", &SolutionBuilderAndSolverType::BuildAndSolve)
          .def("BuildRHSAndSolve", &SolutionBuilderAndSolverType::BuildRHSAndSolve)
          .def("SetUpDofSet", &SolutionBuilderAndSolverType::SetUpDofSet)
          .def("GetDofSet", &SolutionBuilderAndSolverType::GetDofSet, return_internal_reference<>())
          .def("SetUpSystem", &SolutionBuilderAndSolverType::SetUpSystem)
          .def("SetUpSystemMatrices", &SolutionBuilderAndSolverType::SetUpSystemMatrices)
          .def("InitializeSolutionStep", &SolutionBuilderAndSolverType::InitializeSolutionStep)
          .def("FinalizeSolutionStep", &SolutionBuilderAndSolverType::FinalizeSolutionStep)
          .def("CalculateReactions", &SolutionBuilderAndSolverType::CalculateReactions)
          .def("GetEquationSystemSize", &SolutionBuilderAndSolverType::GetEquationSystemSize)
          .def("Clear", &SolutionBuilderAndSolverType::Clear)
          .def("Check", &SolutionBuilderAndSolverType::Check)
          .def("SetEchoLevel", &SolutionBuilderAndSolverType::SetEchoLevel)
          .def("GetEchoLevel", &SolutionBuilderAndSolverType::GetEchoLevel)
          ;

      class_< ReductionBuilderAndSolverType, bases<SolutionBuilderAndSolverType>, boost::noncopyable >
        (
         "ReductionBuilderAndSolver", init< LinearSolverType::Pointer > ()
         );
            
      class_< BlockBuilderAndSolverType, bases<SolutionBuilderAndSolverType>, boost::noncopyable >
        (
         "BlockBuilderAndSolver", init< LinearSolverType::Pointer > ()
         );

      class_< ExplicitBuilderAndSolverType, bases<SolutionBuilderAndSolverType>, boost::noncopyable > 
	(
	 "ExplicitBuilderAndSolver", init<> ()
	 );
      
 
      // class_< ExplicitHamiltonBuilderAndSolverType, bases<BuilderAndSolverType>, boost::noncopyable > 
      //   (
      //    "ExplicitHamiltonBuilderAndSolver", init< LinearSolverType::Pointer > ()
      //    );

      
      //*************************SHCHEME CLASSES****************************

      // Solid Mechanics Base Scheme      
      class_< SolutionSchemeType, bases<Flags>, boost::noncopyable >
          ("SolutionScheme", init< >())
          .def(init < Flags& >())
          .def(init < TimeIntegrationMethodPointerType, Flags& >())
          .def("Initialize", &SolutionSchemeType::Initialize)
          .def("InitializeSolutionStep", &SolutionSchemeType::InitializeSolutionStep)
          .def("FinalizeSolutionStep", &SolutionSchemeType::FinalizeSolutionStep)
          .def("Predict", &SolutionSchemeType::Predict)
          .def("Update", &SolutionSchemeType::Update)
          .def("MoveMesh", &SolutionSchemeType::MoveMesh)
          .def("Check", &SolutionSchemeType::Check)
          ;
      
      // Explicit scheme: Central differences 
      class_< ExplicitCentralDifferencesSchemeType,
	      bases< SolutionSchemeType >,  boost::noncopyable >
	(
          "ExplicitCentralDifferencesScheme", init< Flags& ,const double, const double, const double>() )
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
          .def(init < Flags& >())
      	;
	    
      // Displacement Rotation Static Scheme Type
      class_< DisplacementRotationStaticSchemeType,
      	      bases< SolutionSchemeType >,  boost::noncopyable >
      	(
      	 "DisplacementRotationStaticScheme", init<>() )
          .def(init < Flags& >())
      	;

      // Displacement Newmark Scheme Type
      class_< DisplacementNewmarkSchemeType,
      	      bases< SolutionSchemeType >,  boost::noncopyable >
      	(
      	 "DisplacementNewmarkScheme", init<>() )
          .def(init < Flags& >())
      	;
	    
      // Displacement Rotation Newmark Scheme Type
      class_< DisplacementRotationNewmarkSchemeType,
      	      bases< SolutionSchemeType >,  boost::noncopyable >
      	(
      	 "DisplacementRotationNewmarkScheme", init<>() )
          .def(init < Flags& >())
      	;
      
      // Displacement Bossak Scheme Type
      class_< DisplacementBossakSchemeType,
      	      bases< SolutionSchemeType >,  boost::noncopyable >
      	(
      	 "DisplacementBossakScheme", init<>() )
          .def(init < Flags& >())
      	;
      
      // Displacement Rotation Bossak Scheme Type
      class_< DisplacementRotationBossakSchemeType,
      	      bases< SolutionSchemeType >,  boost::noncopyable >
      	(
      	 "DisplacementRotationBossakScheme", init<>() )
          .def(init < Flags& >())
      	;

      // Displacement Simo Scheme Type
      class_< DisplacementSimoSchemeType,
      	      bases< SolutionSchemeType >,  boost::noncopyable >
      	(
      	 "DisplacementSimoScheme", init<>() )
          .def(init < Flags& >())
       	;

      // Displacement Rotation Simo Scheme Type
      class_< DisplacementRotationSimoSchemeType,
      	      bases< SolutionSchemeType >,  boost::noncopyable >
      	(
      	 "DisplacementRotationSimoScheme", init<>() )
          .def(init < Flags& >())
      	;
      
      // Displacement Rotation Emc Scheme Type
      class_< DisplacementRotationEmcSchemeType,
      	      bases< SolutionSchemeType >,  boost::noncopyable >
      	(
      	 "DisplacementRotationEmcScheme", init<>() )
          .def(init < Flags& >())
      	;


      // Eigensolver Scheme Type
      class_< EigensolverDynamicSchemeType,
	      bases< SolutionSchemeType >, boost::noncopyable >
	(
	 "EigensolverDynamicScheme", init<>() )
          .def(init < Flags& >())
	;

      //*******************CONVERGENCE CRITERIA CLASSES*********************

      // Displacement Convergence Criterion
      class_< DisplacementConvergenceCriterionType,
	      bases< ConvergenceCriteriaType >, boost::noncopyable >
	(
	 "DisplacementConvergenceCriterion", 
	 init<double, double >())
	.def("SetEchoLevel", &DisplacementConvergenceCriterionType::SetEchoLevel)
	;


      //*******************TIME INTEGRATION METHODS*************************

      // Variable type for schemes
      
      //Time Integration Method for schemes
      class_< TimeIntegrationMethodType, TimeIntegrationMethodPointerType, boost::noncopyable >
      	(
      	 "TimeIntegrationMethod", init<>())
	.def("Clone", &TimeIntegrationMethodType::Clone)
	.def("SetVariable", &TimeIntegrationMethodType::SetVariable)
	.def("SetFirstDerivative", &TimeIntegrationMethodType::SetFirstDerivative)
	.def("SetSecondDerivative", &TimeIntegrationMethodType::SetSecondDerivative)
      	.def("SetVariables", &TimeIntegrationMethodType::SetVariables)
      	.def("SetInputVariable", &TimeIntegrationMethodType::SetInputVariable)
	.def("HasStepVariable", &TimeIntegrationMethodType::HasStepVariable)
	.def("SetStepVariable", &TimeIntegrationMethodType::SetStepVariable)
        .def("CalculateParameters", &TimeIntegrationMethodType::CalculateParameters)
      	.def("SetParameters", &TimeIntegrationMethodType::SetParameters)
      	.def("Predict", &TimeIntegrationMethodType::Predict)
	.def(self_ns::str(self))
      	;
      
      // Variable component type for variables integration
      class_< TimeIntegrationContainerType, TimeIntegrationContainerPointerType, boost::noncopyable >
      	(
      	 "TimeIntegrationMethodsContainer", init<>())
        .def("Set", &TimeIntegrationContainerType::Set)
        .def("Get", &TimeIntegrationContainerType::Get)
        .def("Has", &TimeIntegrationContainerType::Has)
        .def(self_ns::str(self))
        DECLARE_HAS_THIS_TYPE_PROCESS_INFO_PYTHON_AS_POINTER(TimeIntegrationContainerType)
        DECLARE_ADD_THIS_TYPE_TO_PROCESS_INFO_PYTHON_AS_POINTER(TimeIntegrationContainerType)
        DECLARE_GET_THIS_TYPE_FROM_PROCESS_INFO_PYTHON_AS_POINTER(TimeIntegrationContainerType)
        ;

      //to define it as a variable 
      class_<Variable<TimeIntegrationContainerPointerType> , bases<VariableData>, boost::noncopyable >("TimeIntegrationMethodsContainerVariable", no_init)
	;
      
      //Time Integration Method for variables
      class_< TimeIntegrationComponentMethodType, TimeIntegrationComponentMethodType::Pointer, boost::noncopyable >
      	(
      	 "TimeIntegrationComponentMethod", init<>())
	.def("Clone", &TimeIntegrationComponentMethodType::Clone)
	.def("SetVariable", &TimeIntegrationComponentMethodType::SetVariable)
	.def("SetFirstDerivative", &TimeIntegrationComponentMethodType::SetFirstDerivative)
	.def("SetSecondDerivative", &TimeIntegrationComponentMethodType::SetSecondDerivative)
      	.def("SetVariables", &TimeIntegrationComponentMethodType::SetVariables)
      	.def("SetInputVariable", &TimeIntegrationComponentMethodType::SetInputVariable)
	.def("HasStepVariable", &TimeIntegrationComponentMethodType::HasStepVariable)
	.def("SetStepVariable", &TimeIntegrationComponentMethodType::SetStepVariable)
        .def("CalculateParameters", &TimeIntegrationComponentMethodType::CalculateParameters)
      	.def("SetParameters", &TimeIntegrationComponentMethodType::SetParameters)
      	.def("Predict", &TimeIntegrationComponentMethodType::Predict)
	.def(self_ns::str(self))
      	;


      class_< StaticMethodType, StaticMethodType::Pointer,
      	      bases< TimeIntegrationComponentMethodType >, boost::noncopyable >
      	(
      	 "StaticMethod", init<>())
       	;
      
      class_< NewmarkMethodType, NewmarkMethodType::Pointer,
      	      bases< TimeIntegrationComponentMethodType >, boost::noncopyable >
      	(
      	 "NewmarkMethod", init<>())
       	;

      class_< BossakMethodType, BossakMethodType::Pointer,
      	      bases< TimeIntegrationComponentMethodType >, boost::noncopyable >
      	(
      	 "BossakMethod", init<>())
      	;
      
      class_< SimoMethodType, SimoMethodType::Pointer,
      	      bases< TimeIntegrationComponentMethodType >, boost::noncopyable >
      	(
      	 "SimoMethod", init<>())
      	;
      
      class_< StaticStepMethodType, StaticStepMethodType::Pointer,
      	      bases< TimeIntegrationComponentMethodType >, boost::noncopyable >
      	(
      	 "StaticStepMethod", init<>())
     	;

      class_< NewmarkStepMethodType, NewmarkStepMethodType::Pointer,
      	      bases< TimeIntegrationComponentMethodType >, boost::noncopyable >
      	(
      	 "NewmarkStepMethod", init<>())
      	;

      class_< BossakStepMethodType, BossakStepMethodType::Pointer,
      	      bases< TimeIntegrationComponentMethodType >, boost::noncopyable >
      	(
      	 "BossakStepMethod", init<>())
      	;

      class_< SimoStepMethodType, SimoStepMethodType::Pointer,
      	      bases< TimeIntegrationComponentMethodType >, boost::noncopyable >
      	(
      	 "SimoStepMethod", init<>())
      	;
      
       class_< EmcStepMethodType, EmcStepMethodType::Pointer,
      	      bases< TimeIntegrationComponentMethodType >, boost::noncopyable >
      	(
      	 "EmcStepMethod", init<>())
      	;

       class_< StaticStepRotationMethodType, StaticStepRotationMethodType::Pointer,
      	      bases< TimeIntegrationComponentMethodType >, boost::noncopyable >
      	(
      	 "StaticStepRotationMethod", init<>())
      	;

       class_< NewmarkStepRotationMethodType, NewmarkStepRotationMethodType::Pointer,
      	      bases< TimeIntegrationComponentMethodType >, boost::noncopyable >
      	(
      	 "NewmarkStepRotationMethod", init<>())
      	;

       class_< BossakStepRotationMethodType, BossakStepRotationMethodType::Pointer,
      	      bases< TimeIntegrationComponentMethodType >, boost::noncopyable >
      	(
      	 "BossakStepRotationMethod", init<>())
      	;

       class_< SimoStepRotationMethodType, SimoStepRotationMethodType::Pointer,
      	      bases< TimeIntegrationComponentMethodType >, boost::noncopyable >
      	(
      	 "SimoStepRotationMethod", init<>())
       	;

       class_< EmcStepRotationMethodType, EmcStepRotationMethodType::Pointer,
      	      bases< TimeIntegrationComponentMethodType >, boost::noncopyable >
      	(
      	 "EmcStepRotationMethod", init<>())
      	;
    }

  }  // namespace Python.

} // Namespace Kratos

