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

// Project includes
#include "spaces/ublas_space.h"
#include "utilities/openmp_utils.h"
#include "custom_python/add_custom_strategies_to_python.h"

// Solution strategies
#include "custom_solvers/solution_strategies/newton_raphson_strategy.hpp"
#include "custom_solvers/solution_strategies/line_search_strategy.hpp"
#include "custom_solvers/solution_strategies/explicit_strategy.hpp" 
#include "custom_solvers/solution_strategies/eigensolver_strategy.hpp"
#include "custom_solvers/solution_strategies/segregated_strategy.hpp"

// to update
#include "custom_solvers/solution_strategies/explicit_hamilton_strategy.hpp"

// Solution builders and solvers
#include "custom_solvers/solution_builders_and_solvers/block_builder_and_solver.hpp"
#include "custom_solvers/solution_builders_and_solvers/reduction_builder_and_solver.hpp"
#include "custom_solvers/solution_builders_and_solvers/explicit_builder_and_solver.hpp"

// to update
#include "custom_solvers/solution_builders_and_solvers/explicit_hamilton_builder_and_solver.hpp"

// Convergence criteria
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "custom_solvers/convergence_criteria/displacement_convergence_criterion.hpp"

// Solution schemes
#include "custom_solvers/solution_schemes/displacement_bossak_scheme.hpp"
#include "custom_solvers/solution_schemes/eigensolver_dynamic_scheme.hpp" 
 
#include "custom_solvers/solution_schemes/explicit_central_differences_scheme.hpp" 
#include "custom_solvers/solution_schemes/explicit_hamilton_scheme.hpp"

#include "custom_solvers/solution_schemes/displacement_rotation_static_scheme.hpp"
#include "custom_solvers/solution_schemes/displacement_rotation_emc_scheme.hpp"
#include "custom_solvers/solution_schemes/displacement_simo_scheme.hpp"
#include "custom_solvers/solution_schemes/displacement_backward_euler_scheme.hpp"
#include "custom_solvers/solution_schemes/displacement_bdf_scheme.hpp"

// Linear solvers
#include "linear_solvers/linear_solver.h"

// Time integration methods
#include "custom_solvers/time_integration_methods/time_integration_methods_container.hpp"


namespace Kratos
{

namespace Python
{
using namespace pybind11;

//base types
typedef UblasSpace<double, CompressedMatrix, Vector>                                               SparseSpaceType;
typedef UblasSpace<double, Matrix, Vector>                                                          LocalSpaceType;
typedef LinearSolver<SparseSpaceType, LocalSpaceType>                                             LinearSolverType;
typedef SolutionStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>                   SolutionStrategyType;
typedef SolutionBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>   SolutionBuilderAndSolverType;
typedef SolutionScheme<SparseSpaceType, LocalSpaceType>                                         SolutionSchemeType;
typedef ConvergenceCriteria<SparseSpaceType, LocalSpaceType>                               ConvergenceCriteriaType; 
typedef SolutionStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>                   SolutionStrategyType;
  
typedef SolutionStrategyType::Pointer                                                      SolutionStrategyPointer;
typedef std::vector<SolutionStrategyType::Pointer>                                     SolutionStrategiesContainer;
  
void Push_Back_Solution_Strategies( SolutionStrategiesContainer& ThisSolutionStrategyContainer,
                                    SolutionStrategyPointer ThisSolutionStrategy )
{ 
  ThisSolutionStrategyContainer.push_back( ThisSolutionStrategy );
} 

void  AddCustomStrategiesToPython(pybind11::module& m)
{

  // Solution strategy types
  typedef SegregatedStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>             SegregatedStrategyType;
  typedef LinearStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>                     LinearStrategyType;
  typedef NewtonRaphsonStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>       NewtonRaphsonStrategyType;
  typedef LineSearchSolutionStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>     LineSearchStrategyType;
  typedef ExplicitSolutionStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>         ExplicitStrategyType;
  typedef ExplicitHamiltonStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> ExplicitHamiltonStrategyType;
  typedef EigensolverStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>           EigensolverStrategyType;

  // Solution builder_and_solver types
  typedef ReductionBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>               ReductionBuilderAndSolverType;
  typedef BlockBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>                       BlockBuilderAndSolverType;
  typedef ExplicitBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>                 ExplicitBuilderAndSolverType;
  typedef ExplicitHamiltonBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> ExplicitHamiltonBuilderAndSolverType;

  // Solution scheme types
  typedef ExplicitCentralDifferencesScheme<SparseSpaceType, LocalSpaceType>    ExplicitCentralDifferencesSchemeType;
  typedef ExplicitHamiltonScheme<SparseSpaceType, LocalSpaceType>                        ExplicitHamiltonSchemeType;
  typedef EigensolverDynamicScheme<SparseSpaceType, LocalSpaceType>                    EigensolverDynamicSchemeType;

  typedef DisplacementStaticScheme<SparseSpaceType, LocalSpaceType>                    DisplacementStaticSchemeType;     
  typedef DisplacementRotationStaticScheme<SparseSpaceType, LocalSpaceType>    DisplacementRotationStaticSchemeType;
      
  typedef DisplacementNewmarkScheme<SparseSpaceType, LocalSpaceType>                  DisplacementNewmarkSchemeType;     
  typedef DisplacementRotationNewmarkScheme<SparseSpaceType, LocalSpaceType>  DisplacementRotationNewmarkSchemeType;
      
  typedef DisplacementBossakScheme<SparseSpaceType, LocalSpaceType>                    DisplacementBossakSchemeType;
  typedef DisplacementRotationBossakScheme<SparseSpaceType, LocalSpaceType>    DisplacementRotationBossakSchemeType;

  typedef DisplacementSimoScheme<SparseSpaceType, LocalSpaceType>                        DisplacementSimoSchemeType;
  typedef DisplacementBackwardEulerScheme<SparseSpaceType, LocalSpaceType>      DisplacementBackwardEulerSchemeType;
  typedef DisplacementBdfScheme<SparseSpaceType, LocalSpaceType>                          DisplacementBdfSchemeType;

  typedef DisplacementRotationSimoScheme<SparseSpaceType, LocalSpaceType>        DisplacementRotationSimoSchemeType;
  typedef DisplacementRotationEmcScheme<SparseSpaceType, LocalSpaceType>          DisplacementRotationEmcSchemeType;
      
  // Custom convergence criterion types
  typedef DisplacementConvergenceCriterion<SparseSpaceType,  LocalSpaceType>   DisplacementConvergenceCriterionType;

  // Time integration methods by variables
  typedef array_1d<double, 3>                                                            VectorType;
  typedef Variable<VectorType>                                                         VariableType; 
  typedef TimeIntegrationMethod<VariableType, VectorType>                 TimeIntegrationMethodType;
  typedef typename TimeIntegrationMethodType::Pointer              TimeIntegrationMethodPointerType;

  // Time integration methods by components
  typedef TimeIntegrationMethodsContainer                              TimeIntegrationContainerType;      
  typedef TimeIntegrationContainerType::Pointer                 TimeIntegrationContainerPointerType;
      
  typedef VariableComponent<VectorComponentAdaptor<VectorType>>               VariableComponentType;
  typedef TimeIntegrationMethod<VariableComponentType, double>   TimeIntegrationComponentMethodType;      
  typedef StaticMethod<VariableComponentType, double>                              StaticMethodType;
  typedef NewmarkMethod<VariableComponentType, double>                            NewmarkMethodType;
  typedef BossakMethod<VariableComponentType, double>                              BossakMethodType;
  typedef SimoMethod<VariableComponentType, double>                                  SimoMethodType;
  typedef BackwardEulerMethod<VariableComponentType, double>                BackwardEulerMethodType;
  typedef BdfMethod<VariableComponentType, double>                                    BdfMethodType;
  
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
  class_<SolverLocalFlags>(m,"SolverLocalFlags")
      .def(init<>())
      .def_readonly_static("INITIALIZED", &SolverLocalFlags::INITIALIZED)
      .def_readonly_static("CONVERGED", &SolverLocalFlags::CONVERGED)
      .def_readonly_static("REFORM_DOFS", &SolverLocalFlags::REFORM_DOFS)
      .def_readonly_static("COMPUTE_REACTIONS", &SolverLocalFlags::COMPUTE_REACTIONS)
      .def_readonly_static("CONSTANT_SYSTEM_MATRIX", &SolverLocalFlags::CONSTANT_SYSTEM_MATRIX)
      .def_readonly_static("RAYLEIGH_DAMPING", &SolverLocalFlags::RAYLEIGH_DAMPING)
      .def_readonly_static("IMPLEX", &SolverLocalFlags::IMPLEX)
      ;

      
  //*************************STRATEGY CLASSES***************************

  // Solid Mechanics Solution Strategies Container
  class_<SolutionStrategiesContainer>(m,"SolutionStragetiesContainer")
      .def(init<>())
      .def("PushBack", Push_Back_Solution_Strategies)
      ;
      
  // Solid Mechanics Base Solution Strategy
  class_<SolutionStrategyType, typename SolutionStrategyType::Pointer, Flags>(m,"SolutionStrategy")
      .def(init<ModelPart&>())
      .def(init<ModelPart&, Flags&>())
      .def("InitializeSolutionStep", &SolutionStrategyType::InitializeSolutionStep)
      .def("FinalizeSolutionStep", &SolutionStrategyType::FinalizeSolutionStep)
      .def("SolveSolutionStep", &SolutionStrategyType::SolveSolutionStep )
      .def("Solve", &SolutionStrategyType::Solve)
      .def("Check", &SolutionStrategyType::Check)
      .def("Clear", &SolutionStrategyType::Clear)
      .def("SetOptions", &SolutionStrategyType::SetOptions)
      .def("GetOptions", &SolutionStrategyType::GetOptions, return_value_policy::reference_internal)
      .def("SetEchoLevel", &SolutionStrategyType::SetEchoLevel)
      .def("GetEchoLevel", &SolutionStrategyType::GetEchoLevel)
      ;

  // Solid Mechanics Segregated Strategy
  class_<SegregatedStrategyType, typename SegregatedStrategyType::Pointer, SolutionStrategyType>(m,"SegregatedStrategy")
      .def(init<ModelPart&, Flags&>())
      .def(init<ModelPart&, Flags&,  SolutionStrategiesContainer&>())
      .def("AddStrategy", &SegregatedStrategyType::AddStrategy)
      .def("GetStrategy", &SegregatedStrategyType::GetStrategy)
      ;
      
  // Solid Mechanics Linear Strategy
  class_<LinearStrategyType, typename LinearStrategyType::Pointer, SolutionStrategyType>(m,"LinearStrategy")
      .def(init<ModelPart&, SolutionSchemeType::Pointer, SolutionBuilderAndSolverType::Pointer, Flags&>())
      .def(init<ModelPart&, SolutionSchemeType::Pointer, LinearSolverType::Pointer, Flags&>())
      ;
      
  // Solid Mechanics Newton Raphson Strategy
  class_<NewtonRaphsonStrategyType, typename NewtonRaphsonStrategyType::Pointer, LinearStrategyType>(m,"NewtonRaphsonStrategy")
      .def(init<ModelPart&, SolutionSchemeType::Pointer, SolutionBuilderAndSolverType::Pointer, ConvergenceCriteriaType::Pointer, Flags&, unsigned int>())
      .def(init<ModelPart&, SolutionSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, Flags&, unsigned int>())
      .def("SetMaxIterationNumber", &NewtonRaphsonStrategyType::SetMaxIterationNumber)
      .def("GetMaxIterationNumber", &NewtonRaphsonStrategyType::GetMaxIterationNumber)  
      ;

  // Solid Mechanics Newton Raphson Line Search Strategy
  class_<LineSearchStrategyType, typename LineSearchStrategyType::Pointer, NewtonRaphsonStrategyType>(m,"LineSearchStrategy")
      .def(init<ModelPart&, SolutionSchemeType::Pointer, SolutionBuilderAndSolverType::Pointer, ConvergenceCriteriaType::Pointer, Flags&, unsigned int>())
      .def(init<ModelPart&, SolutionSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, Flags&, unsigned int>())
      ;
      
  // Solid Mechanics Explicit Strategy
  class_<ExplicitStrategyType, typename ExplicitStrategyType::Pointer, SolutionStrategyType>(m,"ExplicitStrategy")
      .def(init<ModelPart&, SolutionSchemeType::Pointer, Flags&>())
      ;
         

  // Eigensolver Strategy
  class_<EigensolverStrategyType, typename EigensolverStrategyType::Pointer, SolutionStrategyType>(m,"EigensolverStrategy")
      .def(init<ModelPart&, SolutionSchemeType::Pointer, SolutionBuilderAndSolverType::Pointer, Flags&, bool>())
      ;

      
  // // Explicit Hamilton Estrategy for Explicit Beam solution
  // class_<ExplicitHamiltonStrategyType, BaseSolvingStrategyType>(m,"ExplicitHamiltonStrategy")
  //   .def(init<ModelPart&, SolutionSchemeType::Pointer, LinearSolverType::Pointer, bool, bool, bool>())	
  //   ;

      
  //*******************BUILDER AND SOLVER CLASSES***********************

  // Solid Mechanics Base Builder and Solver
  class_<SolutionBuilderAndSolverType, typename SolutionBuilderAndSolverType::Pointer, Flags>
      (m,"SolutionBuilderAndSolver")
      .def(init<LinearSolverType::Pointer> ())
      .def(init<>())
      .def("BuildLHS", &SolutionBuilderAndSolverType::BuildLHS)
      .def("BuildRHS", &SolutionBuilderAndSolverType::BuildRHS)
      .def("Build", &SolutionBuilderAndSolverType::Build)
      .def("SystemSolve", &SolutionBuilderAndSolverType::SystemSolve)
      .def("BuildAndSolve", &SolutionBuilderAndSolverType::BuildAndSolve)
      .def("BuildRHSAndSolve", &SolutionBuilderAndSolverType::BuildRHSAndSolve)
      .def("SetUpDofSet", &SolutionBuilderAndSolverType::SetUpDofSet)
      .def("GetDofSet", &SolutionBuilderAndSolverType::GetDofSet, return_value_policy::reference_internal)
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

  class_<ReductionBuilderAndSolverType, typename ReductionBuilderAndSolverType::Pointer, SolutionBuilderAndSolverType>(m,"ReductionBuilderAndSolver")
      .def(init<LinearSolverType::Pointer>())
      ;
            
  class_<BlockBuilderAndSolverType, typename BlockBuilderAndSolverType::Pointer, SolutionBuilderAndSolverType>(m,"BlockBuilderAndSolver")
      .def(init<LinearSolverType::Pointer>())
      ;

  class_<ExplicitBuilderAndSolverType, typename ExplicitBuilderAndSolverType::Pointer, SolutionBuilderAndSolverType>(m,"ExplicitBuilderAndSolver")
      .def(init<>())
      ;
      
 
  // class_<ExplicitHamiltonBuilderAndSolverType, typename ExplicitHamiltonBuilderAndSolverType::Pointer, BuilderAndSolverType>(m,"ExplicitHamiltonBuilderAndSolver")
  //   .def(init<LinearSolverType::Pointer> ())
  //;

      
  //*************************SHCHEME CLASSES****************************

  // Solid Mechanics Base Scheme      
  class_<SolutionSchemeType, typename SolutionSchemeType::Pointer, Flags>(m,"SolutionScheme")
      .def(init<>())
      .def(init<Flags&>())
      .def(init<TimeIntegrationMethodPointerType, Flags&>())
      .def("Initialize", &SolutionSchemeType::Initialize)
      .def("InitializeSolutionStep", &SolutionSchemeType::InitializeSolutionStep)
      .def("FinalizeSolutionStep", &SolutionSchemeType::FinalizeSolutionStep)
      .def("Update", &SolutionSchemeType::Update)
      .def("MoveMesh", &SolutionSchemeType::MoveMesh)
      .def("Check", &SolutionSchemeType::Check)
      ;
      
  // Explicit scheme: Central differences 
  class_<ExplicitCentralDifferencesSchemeType, typename ExplicitCentralDifferencesSchemeType::Pointer, SolutionSchemeType>(m,"ExplicitCentralDifferencesScheme")
      .def(init<Flags& ,const double, const double, const double>())
      ;

  // Displacement Static Scheme Type
  class_<DisplacementStaticSchemeType, typename DisplacementStaticSchemeType::Pointer, SolutionSchemeType>(m,"DisplacementStaticScheme")
      .def(init<>())
      .def(init<Flags&>())
      ;
	    
  // Displacement Rotation Static Scheme Type
  class_<DisplacementRotationStaticSchemeType, typename DisplacementRotationStaticSchemeType::Pointer, SolutionSchemeType>(m,"DisplacementRotationStaticScheme")
      .def(init<>())
      .def(init<Flags&>())
      ;

  // Displacement Newmark Scheme Type
  class_<DisplacementNewmarkSchemeType, typename DisplacementNewmarkSchemeType::Pointer, SolutionSchemeType>(m,"DisplacementNewmarkScheme")
      .def(init<>())
      .def(init<Flags&>())
      ;
	    
  // Displacement Rotation Newmark Scheme Type
  class_<DisplacementRotationNewmarkSchemeType,  typename DisplacementRotationNewmarkSchemeType::Pointer, SolutionSchemeType>(m,"DisplacementRotationNewmarkScheme")
      .def(init<>())
      .def(init<Flags&>())
      ;
      
  // Displacement Bossak Scheme Type
  class_<DisplacementBossakSchemeType, typename DisplacementBossakSchemeType::Pointer, SolutionSchemeType>(m,"DisplacementBossakScheme")
      .def(init<>())
      .def(init<Flags&>())
      ;
      
  // Displacement Rotation Bossak Scheme Type
  class_<DisplacementRotationBossakSchemeType, typename DisplacementRotationBossakSchemeType::Pointer, SolutionSchemeType>(m,"DisplacementRotationBossakScheme")
      .def(init<>())
      .def(init<Flags&>())
      ;

  // Displacement Simo Scheme Type
  class_<DisplacementSimoSchemeType,  typename DisplacementSimoSchemeType::Pointer, SolutionSchemeType>(m,"DisplacementSimoScheme")
      .def(init<>())
      .def(init<Flags&>())
      ;

  // Displacement Backward Euler Scheme Type
  class_<DisplacementBackwardEulerSchemeType,  typename DisplacementBackwardEulerSchemeType::Pointer, SolutionSchemeType>(m,"DisplacementBackwardEulerScheme")
      .def(init<>())
      .def(init<Flags&>())
      ;

  // Displacement Bdf Scheme Type
  class_<DisplacementBdfSchemeType,  typename DisplacementBdfSchemeType::Pointer, SolutionSchemeType>(m,"DisplacementBdfScheme")
      .def(init<>())
      .def(init<Flags&>())
      ;
  
  // Displacement Rotation Simo Scheme Type
  class_<DisplacementRotationSimoSchemeType, typename DisplacementRotationSimoSchemeType::Pointer, SolutionSchemeType>(m,"DisplacementRotationSimoScheme")
      .def(init<>())
      .def(init<Flags&>())
      ;
      
  // Displacement Rotation Emc Scheme Type
  class_<DisplacementRotationEmcSchemeType,  typename DisplacementRotationEmcSchemeType::Pointer, SolutionSchemeType>(m,"DisplacementRotationEmcScheme")
      .def(init<>())
      .def(init<Flags&>())
      ;


  // Eigensolver Scheme Type
  class_<EigensolverDynamicSchemeType, typename EigensolverDynamicSchemeType::Pointer, SolutionSchemeType>(m,"EigensolverDynamicScheme")
      .def(init<>())
      .def(init<Flags&>())
      ;


  // // Explicit Hamilton Scheme Type
  // class_<ExplicitHamiltonSchemeType, typename ExplicitHamiltonSchemeType::Pointer, SolutionSchemeType> (m,"ExplicitHamiltonScheme")
  //  .def(init<double, double, double, bool>())	
  //   ;
  
  //*******************CONVERGENCE CRITERIA CLASSES*********************

  // Displacement Convergence Criterion
  class_<DisplacementConvergenceCriterionType, typename DisplacementConvergenceCriterionType::Pointer, ConvergenceCriteriaType>
      (m,"DisplacementConvergenceCriterion")
      .def( init<double, double>())
      .def("SetEchoLevel", &DisplacementConvergenceCriterionType::SetEchoLevel)
      ;


  //*******************TIME INTEGRATION METHODS*************************

  // Variable type for schemes
      
  //Time Integration Method for schemes
  class_<TimeIntegrationMethodType, TimeIntegrationMethodPointerType>(m,"TimeIntegrationMethod")
      .def(init<>())
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
      .def("Assign", &TimeIntegrationMethodType::Assign)
      .def("Predict", &TimeIntegrationMethodType::Predict)
      .def("__repr__", &TimeIntegrationMethodType::Info)
      ;
      
  // Variable component type for variables integration
  class_<TimeIntegrationContainerType, TimeIntegrationContainerPointerType>(m,"TimeIntegrationMethodsContainer")
      .def(init<>())
      .def("Set", &TimeIntegrationContainerType::Set)
      .def("Get", &TimeIntegrationContainerType::Get)
      .def("Has", &TimeIntegrationContainerType::Has)
      .def("__repr__", &TimeIntegrationContainerType::Info)
      DECLARE_HAS_THIS_TYPE_PROCESS_INFO_PYTHON_AS_POINTER(TimeIntegrationContainerType)
      DECLARE_ADD_THIS_TYPE_TO_PROCESS_INFO_PYTHON_AS_POINTER(TimeIntegrationContainerType)
      DECLARE_GET_THIS_TYPE_FROM_PROCESS_INFO_PYTHON_AS_POINTER(TimeIntegrationContainerType)
      ;

  //to define it as a variable 
  class_<Variable<TimeIntegrationContainerPointerType>, VariableData>(m,"TimeIntegrationMethodsContainerVariable")
      ;
      
  //Time Integration Method for variables
  class_<TimeIntegrationComponentMethodType, typename TimeIntegrationComponentMethodType::Pointer>(m,"TimeIntegrationComponentMethod")
      .def(init<>())
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
      .def("Assign", &TimeIntegrationComponentMethodType::Assign)      
      .def("Predict", &TimeIntegrationComponentMethodType::Predict)
      .def("__repr__", &TimeIntegrationComponentMethodType::Info)
      ;


  class_<StaticMethodType, typename StaticMethodType::Pointer,
         TimeIntegrationComponentMethodType>(m,"StaticMethod")
      .def(init<>())
      ;
      
  class_<NewmarkMethodType, typename NewmarkMethodType::Pointer,
         TimeIntegrationComponentMethodType>(m,"NewmarkMethod")
      .def(init<>())
      ;

  class_<BossakMethodType, typename BossakMethodType::Pointer,
         TimeIntegrationComponentMethodType>(m,"BossakMethod")
      .def(init<>())
      ;
      
  class_<SimoMethodType, typename SimoMethodType::Pointer,
         TimeIntegrationComponentMethodType>(m,"SimoMethod")
      .def(init<>())
      ;
      
  class_<BackwardEulerMethodType, typename BackwardEulerMethodType::Pointer,
         TimeIntegrationComponentMethodType>(m,"BackwardEulerMethod")
      .def(init<>())
      ;

  class_<BdfMethodType, typename BdfMethodType::Pointer,
         TimeIntegrationComponentMethodType>(m,"BdfMethod")
      .def(init<>())
      ;
  
  class_<StaticStepMethodType, typename StaticStepMethodType::Pointer,
         TimeIntegrationComponentMethodType>(m,"StaticStepMethod")
      .def(init<>())
      ;

  class_<NewmarkStepMethodType, typename NewmarkStepMethodType::Pointer,
         TimeIntegrationComponentMethodType>(m,"NewmarkStepMethod")
      .def(init<>())
      ;

  class_<BossakStepMethodType, typename BossakStepMethodType::Pointer,
         TimeIntegrationComponentMethodType>(m,"BossakStepMethod")
      .def(init<>())
      ;

  class_<SimoStepMethodType, typename SimoStepMethodType::Pointer,
         TimeIntegrationComponentMethodType>(m,"SimoStepMethod")
      .def(init<>())
      ;
      
  class_<EmcStepMethodType, typename EmcStepMethodType::Pointer,
         TimeIntegrationComponentMethodType>(m,"EmcStepMethod")
      .def(init<>())
      ;

  class_<StaticStepRotationMethodType, typename StaticStepRotationMethodType::Pointer,
         TimeIntegrationComponentMethodType>(m,"StaticStepRotationMethod")
      .def(init<>())
      ;

  class_<NewmarkStepRotationMethodType, typename NewmarkStepRotationMethodType::Pointer,
         TimeIntegrationComponentMethodType>(m,"NewmarkStepRotationMethod")
      .def(init<>())
      ;

  class_<BossakStepRotationMethodType, typename BossakStepRotationMethodType::Pointer,
         TimeIntegrationComponentMethodType>(m,"BossakStepRotationMethod")
      .def(init<>())
      ;

  class_<SimoStepRotationMethodType, typename SimoStepRotationMethodType::Pointer,
         TimeIntegrationComponentMethodType>(m,"SimoStepRotationMethod")
      .def(init<>())
      ;

  class_<EmcStepRotationMethodType, typename EmcStepRotationMethodType::Pointer,
         TimeIntegrationComponentMethodType>(m,"EmcStepRotationMethod")
      .def(init<>())
      ;
}

}  // namespace Python.

} // Namespace Kratos

