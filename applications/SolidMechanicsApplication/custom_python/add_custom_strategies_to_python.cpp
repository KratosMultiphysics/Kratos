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
#include "custom_solvers/solution_schemes/static_scheme.hpp"

#include "custom_solvers/solution_schemes/displacement_bossak_scheme.hpp"
#include "custom_solvers/solution_schemes/eigensolver_dynamic_scheme.hpp"

#include "custom_solvers/solution_schemes/explicit_central_differences_scheme.hpp"
#include "custom_solvers/solution_schemes/explicit_hamilton_scheme.hpp"

#include "custom_solvers/solution_schemes/displacement_rotation_static_scheme.hpp"
#include "custom_solvers/solution_schemes/displacement_rotation_emc_scheme.hpp"
#include "custom_solvers/solution_schemes/displacement_simo_scheme.hpp"
#include "custom_solvers/solution_schemes/displacement_backward_euler_scheme.hpp"

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

  typedef StaticScheme<SparseSpaceType, LocalSpaceType>                                            StaticSchemeType;
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
  typedef array_1d<double, 3>                                                                  VectorType;
  typedef Variable<VectorType>                                                               VariableType;
  typedef TimeIntegrationMethod<VariableType, VectorType>                 TimeIntegrationMethodVectorType;
  
  typedef std::vector<TimeIntegrationMethodVectorType::Pointer>              TimeIntegrationMethodsVector;

  typedef StaticMethod<VariableType, VectorType>                                   StaticMethodVectorType;
  typedef NewmarkMethod<VariableType, VectorType>                                 NewmarkMethodVectorType;
  typedef BossakMethod<VariableType, VectorType>                                   BossakMethodVectorType;
  typedef SimoMethod<VariableType, VectorType>                                       SimoMethodVectorType;
  typedef BackwardEulerMethod<VariableType, VectorType>                     BackwardEulerMethodVectorType;
  typedef BdfMethod<VariableType, VectorType>                                         BdfMethodVectorType;

  typedef StaticStepMethod<VariableType, VectorType>                           StaticStepMethodVectorType;
  typedef NewmarkStepMethod<VariableType, VectorType>                         NewmarkStepMethodVectorType;
  typedef BossakStepMethod<VariableType, VectorType>                           BossakStepMethodVectorType;
  typedef SimoStepMethod<VariableType, VectorType>                               SimoStepMethodVectorType;
  typedef EmcStepMethod<VariableType, VectorType>                                 EmcStepMethodVectorType;

  typedef StaticStepRotationMethod<VariableType, VectorType>           StaticStepRotationMethodVectorType;
  typedef NewmarkStepRotationMethod<VariableType, VectorType>         NewmarkStepRotationMethodVectorType;
  typedef BossakStepRotationMethod<VariableType, VectorType>           BossakStepRotationMethodVectorType;
  typedef SimoStepRotationMethod<VariableType, VectorType>               SimoStepRotationMethodVectorType;
  typedef EmcStepRotationMethod<VariableType, VectorType>                 EmcStepRotationMethodVectorType;
       
  
  // Time integration methods by components 
  typedef VariableComponent<VectorComponentAdaptor<VectorType>>                     VariableComponentType;
  typedef TimeIntegrationMethod<VariableComponentType, double>            TimeIntegrationMethodScalarType;

  typedef TimeIntegrationMethodsContainer                                    TimeIntegrationContainerType;
  typedef TimeIntegrationContainerType::Pointer                       TimeIntegrationContainerPointerType;

  typedef StaticMethod<VariableComponentType, double>                              StaticMethodScalarType;
  typedef NewmarkMethod<VariableComponentType, double>                            NewmarkMethodScalarType;
  typedef BossakMethod<VariableComponentType, double>                              BossakMethodScalarType;
  typedef SimoMethod<VariableComponentType, double>                                  SimoMethodScalarType;
  typedef BackwardEulerMethod<VariableComponentType, double>                BackwardEulerMethodScalarType;
  typedef BdfMethod<VariableComponentType, double>                                    BdfMethodScalarType;

  typedef StaticStepMethod<VariableComponentType, double>                      StaticStepMethodScalarType;
  typedef NewmarkStepMethod<VariableComponentType, double>                    NewmarkStepMethodScalarType;
  typedef BossakStepMethod<VariableComponentType, double>                      BossakStepMethodScalarType;
  typedef SimoStepMethod<VariableComponentType, double>                          SimoStepMethodScalarType;
  typedef EmcStepMethod<VariableComponentType, double>                            EmcStepMethodScalarType;

  typedef StaticStepRotationMethod<VariableComponentType, double>      StaticStepRotationMethodScalarType;
  typedef NewmarkStepRotationMethod<VariableComponentType, double>    NewmarkStepRotationMethodScalarType;
  typedef BossakStepRotationMethod<VariableComponentType, double>      BossakStepRotationMethodScalarType;
  typedef SimoStepRotationMethod<VariableComponentType, double>          SimoStepRotationMethodScalarType;
  typedef EmcStepRotationMethod<VariableComponentType, double>            EmcStepRotationMethodScalarType;
  
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
      .def(init<TimeIntegrationMethodsVector&, Flags&>())
      .def(init<TimeIntegrationMethodsVector&>())
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
  class_<StaticSchemeType, typename StaticSchemeType::Pointer, SolutionSchemeType>(m,"StaticScheme")
      .def(init<TimeIntegrationMethodsVector&, Flags&>())
      .def(init<TimeIntegrationMethodsVector&>())
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
  class_<TimeIntegrationMethodVectorType, TimeIntegrationMethodVectorType::Pointer>(m,"TimeIntegration")
      .def(init<>())
      .def(init<const VariableType&>())
      .def(init<const VariableType&, const VariableType&, const VariableType&>())
      .def(init<const VariableType&, const VariableType&, const VariableType&, const VariableType&>())
      .def("Clone", &TimeIntegrationMethodVectorType::Clone)
      .def("SetVariable", &TimeIntegrationMethodVectorType::SetVariable)
      .def("SetFirstDerivative", &TimeIntegrationMethodVectorType::SetFirstDerivative)
      .def("SetSecondDerivative", &TimeIntegrationMethodVectorType::SetSecondDerivative)
      .def("SetVariables", &TimeIntegrationMethodVectorType::SetVariables)
      .def("SetInputVariable", &TimeIntegrationMethodVectorType::SetInputVariable)
      .def("HasStepVariable", &TimeIntegrationMethodVectorType::HasStepVariable)
      .def("SetStepVariable", &TimeIntegrationMethodVectorType::SetStepVariable)
      .def("CalculateParameters", &TimeIntegrationMethodVectorType::CalculateParameters)
      .def("SetParameters", &TimeIntegrationMethodVectorType::SetParameters)
      .def("Assign", &TimeIntegrationMethodVectorType::Assign)
      .def("Predict", &TimeIntegrationMethodVectorType::Predict)
      .def("__repr__", &TimeIntegrationMethodVectorType::Info)
      ;


  class_<StaticMethodVectorType, typename StaticMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"StaticIntegration")
      .def(init<const VariableType&>())
      ;

  class_<NewmarkMethodVectorType, typename NewmarkMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"NewmarkIntegration")
      .def(init<const VariableType&, const VariableType&, const VariableType&>())
      .def(init<const VariableType&, const VariableType&, const VariableType&, const VariableType&>())
      ;

  class_<BossakMethodVectorType, typename BossakMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"BossakIntegration")
      .def(init<const VariableType&, const VariableType&, const VariableType&>())
      .def(init<const VariableType&, const VariableType&, const VariableType&, const VariableType&>())
      ;

  class_<SimoMethodVectorType, typename SimoMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"SimoIntegration")
      .def(init<const VariableType&, const VariableType&, const VariableType&>())
      .def(init<const VariableType&, const VariableType&, const VariableType&, const VariableType&>())
      ;

  class_<BackwardEulerMethodVectorType, typename BackwardEulerMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"BackwardEulerIntegration")
      .def(init<const VariableType&, const VariableType&, const VariableType&>())
      .def(init<const VariableType&, const VariableType&, const VariableType&, const VariableType&>())
      ;

  class_<BdfMethodVectorType, typename BdfMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"BdfIntegration")
      .def(init<const VariableType&, const VariableType&, const VariableType&>())
      .def(init<const VariableType&, const VariableType&, const VariableType&, const VariableType&>())
      ;

  class_<StaticStepMethodVectorType, typename StaticStepMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"StaticStepIntegration")
      .def(init<const VariableType&>())
      ;

  class_<NewmarkStepMethodVectorType, typename NewmarkStepMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"NewmarkStepIntegration")
      .def(init<const VariableType&, const VariableType&, const VariableType&>())
      .def(init<const VariableType&, const VariableType&, const VariableType&, const VariableType&>())
      ;

  class_<BossakStepMethodVectorType, typename BossakStepMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"BossakStepIntegration")
      .def(init<const VariableType&, const VariableType&, const VariableType&>())
      .def(init<const VariableType&, const VariableType&, const VariableType&, const VariableType&>())
      ;

  class_<SimoStepMethodVectorType, typename SimoStepMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"SimoStepIntegration")
      .def(init<const VariableType&, const VariableType&, const VariableType&>())
      .def(init<const VariableType&, const VariableType&, const VariableType&, const VariableType&>())
      ;

  class_<EmcStepMethodVectorType, typename EmcStepMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"EmcStepIntegration")
      .def(init<const VariableType&, const VariableType&, const VariableType&>())
      .def(init<const VariableType&, const VariableType&, const VariableType&, const VariableType&>())
      ;

  class_<StaticStepRotationMethodVectorType, typename StaticStepRotationMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"StaticStepRotationIntegration")
      .def(init<const VariableType&>())
      ;

  class_<NewmarkStepRotationMethodVectorType, typename NewmarkStepRotationMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"NewmarkStepRotationIntegration")
      .def(init<const VariableType&, const VariableType&, const VariableType&>())
      .def(init<const VariableType&, const VariableType&, const VariableType&, const VariableType&>())
      ;

  class_<BossakStepRotationMethodVectorType, typename BossakStepRotationMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"BossakStepRotationIntegration")
      .def(init<const VariableType&, const VariableType&, const VariableType&>())
      .def(init<const VariableType&, const VariableType&, const VariableType&, const VariableType&>())
      ;

  class_<SimoStepRotationMethodVectorType, typename SimoStepRotationMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"SimoStepRotationIntegration")
      .def(init<const VariableType&, const VariableType&, const VariableType&>())
      .def(init<const VariableType&, const VariableType&, const VariableType&, const VariableType&>())
      ;

  class_<EmcStepRotationMethodVectorType, typename EmcStepRotationMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"EmcStepRotationIntegration")
      .def(init<const VariableType&, const VariableType&, const VariableType&>())
      .def(init<const VariableType&, const VariableType&, const VariableType&, const VariableType&>())
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
  class_<TimeIntegrationMethodScalarType, typename TimeIntegrationMethodScalarType::Pointer>(m,"TimeIntegrationMethod")
      .def(init<>())
      .def("Clone", &TimeIntegrationMethodScalarType::Clone)
      .def("SetVariable", &TimeIntegrationMethodScalarType::SetVariable)
      .def("SetFirstDerivative", &TimeIntegrationMethodScalarType::SetFirstDerivative)
      .def("SetSecondDerivative", &TimeIntegrationMethodScalarType::SetSecondDerivative)
      .def("SetVariables", &TimeIntegrationMethodScalarType::SetVariables)
      .def("SetInputVariable", &TimeIntegrationMethodScalarType::SetInputVariable)
      .def("HasStepVariable", &TimeIntegrationMethodScalarType::HasStepVariable)
      .def("SetStepVariable", &TimeIntegrationMethodScalarType::SetStepVariable)
      .def("CalculateParameters", &TimeIntegrationMethodScalarType::CalculateParameters)
      .def("SetParameters", &TimeIntegrationMethodScalarType::SetParameters)
      .def("Assign", &TimeIntegrationMethodScalarType::Assign)
      .def("Predict", &TimeIntegrationMethodScalarType::Predict)
      .def("__repr__", &TimeIntegrationMethodScalarType::Info)
      ;


  class_<StaticMethodScalarType, typename StaticMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"StaticMethod")
      .def(init<>())
      ;

  class_<NewmarkMethodScalarType, typename NewmarkMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"NewmarkMethod")
      .def(init<>())
      ;

  class_<BossakMethodScalarType, typename BossakMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"BossakMethod")
      .def(init<>())
      ;

  class_<SimoMethodScalarType, typename SimoMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"SimoMethod")
      .def(init<>())
      ;

  class_<BackwardEulerMethodScalarType, typename BackwardEulerMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"BackwardEulerMethod")
      .def(init<>())
      ;

  class_<BdfMethodScalarType, typename BdfMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"BdfMethod")
      .def(init<>())
      ;

  class_<StaticStepMethodScalarType, typename StaticStepMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"StaticStepMethod")
      .def(init<>())
      ;

  class_<NewmarkStepMethodScalarType, typename NewmarkStepMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"NewmarkStepMethod")
      .def(init<>())
      ;

  class_<BossakStepMethodScalarType, typename BossakStepMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"BossakStepMethod")
      .def(init<>())
      ;

  class_<SimoStepMethodScalarType, typename SimoStepMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"SimoStepMethod")
      .def(init<>())
      ;

  class_<EmcStepMethodScalarType, typename EmcStepMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"EmcStepMethod")
      .def(init<>())
      ;

  class_<StaticStepRotationMethodScalarType, typename StaticStepRotationMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"StaticStepRotationMethod")
      .def(init<>())
      ;

  class_<NewmarkStepRotationMethodScalarType, typename NewmarkStepRotationMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"NewmarkStepRotationMethod")
      .def(init<>())
      ;

  class_<BossakStepRotationMethodScalarType, typename BossakStepRotationMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"BossakStepRotationMethod")
      .def(init<>())
      ;

  class_<SimoStepRotationMethodScalarType, typename SimoStepRotationMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"SimoStepRotationMethod")
      .def(init<>())
      ;

  class_<EmcStepRotationMethodScalarType, typename EmcStepRotationMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"EmcStepRotationMethod")
      .def(init<>())
      ;
}

}  // namespace Python.

} // Namespace Kratos
