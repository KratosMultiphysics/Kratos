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

  // Time integration methods for vectors
  typedef array_1d<double, 3>                                                                      VectorType;
  typedef Variable<VectorType>                                                             VariableVectorType;
  typedef TimeIntegrationMethod<VariableVectorType, VectorType>               TimeIntegrationMethodVectorType;
  
  typedef std::vector<TimeIntegrationMethodVectorType::Pointer>                  TimeIntegrationMethodsVector;
  typedef TimeIntegrationMethodsContainer<VariableVectorType, double>      VectorTimeIntegrationContainerType;
  typedef VectorTimeIntegrationContainerType::Pointer               VectorTimeIntegrationContainerPointerType;
  
  typedef StaticMethod<VariableVectorType, VectorType>                                 StaticMethodVectorType;
  typedef NewmarkMethod<VariableVectorType, VectorType>                               NewmarkMethodVectorType;
  typedef BossakMethod<VariableVectorType, VectorType>                                 BossakMethodVectorType;
  typedef SimoMethod<VariableVectorType, VectorType>                                     SimoMethodVectorType;
  typedef BackwardEulerMethod<VariableVectorType, VectorType>                   BackwardEulerMethodVectorType;
  typedef BdfMethod<VariableVectorType, VectorType>                                       BdfMethodVectorType;

  typedef StaticStepMethod<VariableVectorType, VectorType>                         StaticStepMethodVectorType;
  typedef NewmarkStepMethod<VariableVectorType, VectorType>                       NewmarkStepMethodVectorType;
  typedef BossakStepMethod<VariableVectorType, VectorType>                         BossakStepMethodVectorType;
  typedef SimoStepMethod<VariableVectorType, VectorType>                             SimoStepMethodVectorType;
  typedef EmcStepMethod<VariableVectorType, VectorType>                               EmcStepMethodVectorType;

  typedef StaticStepRotationMethod<VariableVectorType, VectorType>         StaticStepRotationMethodVectorType;
  typedef NewmarkStepRotationMethod<VariableVectorType, VectorType>       NewmarkStepRotationMethodVectorType;
  typedef BossakStepRotationMethod<VariableVectorType, VectorType>         BossakStepRotationMethodVectorType;
  typedef SimoStepRotationMethod<VariableVectorType, VectorType>             SimoStepRotationMethodVectorType;
  typedef EmcStepRotationMethod<VariableVectorType, VectorType>               EmcStepRotationMethodVectorType;
       
  // Time integration methods for scalars 
  typedef Variable<double>                                                                 VariableScalarType;
  typedef TimeIntegrationMethod<VariableScalarType, double>                   TimeIntegrationMethodScalarType;

  typedef TimeIntegrationMethodsContainer<VariableScalarType, double>      ScalarTimeIntegrationContainerType;
  typedef ScalarTimeIntegrationContainerType::Pointer               ScalarTimeIntegrationContainerPointerType;

  typedef StaticMethod<VariableScalarType, double>                                     StaticMethodScalarType;
  typedef NewmarkMethod<VariableScalarType, double>                                   NewmarkMethodScalarType;
  typedef BossakMethod<VariableScalarType, double>                                     BossakMethodScalarType;
  typedef SimoMethod<VariableScalarType, double>                                         SimoMethodScalarType;
  typedef BackwardEulerMethod<VariableScalarType, double>                       BackwardEulerMethodScalarType;
  typedef BdfMethod<VariableScalarType, double>                                           BdfMethodScalarType;

  typedef StaticStepMethod<VariableScalarType, double>                             StaticStepMethodScalarType;
  typedef NewmarkStepMethod<VariableScalarType, double>                           NewmarkStepMethodScalarType;
  typedef BossakStepMethod<VariableScalarType, double>                             BossakStepMethodScalarType;
  typedef SimoStepMethod<VariableScalarType, double>                                 SimoStepMethodScalarType;
  typedef EmcStepMethod<VariableScalarType, double>                                   EmcStepMethodScalarType;

  typedef StaticStepRotationMethod<VariableScalarType, double>             StaticStepRotationMethodScalarType;
  typedef NewmarkStepRotationMethod<VariableScalarType, double>           NewmarkStepRotationMethodScalarType;
  typedef BossakStepRotationMethod<VariableScalarType, double>             BossakStepRotationMethodScalarType;
  typedef SimoStepRotationMethod<VariableScalarType, double>                 SimoStepRotationMethodScalarType;
  typedef EmcStepRotationMethod<VariableScalarType, double>                   EmcStepRotationMethodScalarType;

  
  // Time integration methods for components
  typedef VariableComponent<VectorComponentAdaptor<VectorType>>                        VariableComponentType;
  typedef TimeIntegrationMethod<VariableComponentType, double>            TimeIntegrationMethodComponentType;

  typedef TimeIntegrationMethodsContainer<VariableComponentType, double> ComponentTimeIntegrationContainerType;
  typedef ComponentTimeIntegrationContainerType::Pointer          ComponentTimeIntegrationContainerPointerType;

  typedef StaticMethod<VariableComponentType, double>                              StaticMethodComponentType;
  typedef NewmarkMethod<VariableComponentType, double>                            NewmarkMethodComponentType;
  typedef BossakMethod<VariableComponentType, double>                              BossakMethodComponentType;
  typedef SimoMethod<VariableComponentType, double>                                  SimoMethodComponentType;
  typedef BackwardEulerMethod<VariableComponentType, double>                BackwardEulerMethodComponentType;
  typedef BdfMethod<VariableComponentType, double>                                    BdfMethodComponentType;

  typedef StaticStepMethod<VariableComponentType, double>                      StaticStepMethodComponentType;
  typedef NewmarkStepMethod<VariableComponentType, double>                    NewmarkStepMethodComponentType;
  typedef BossakStepMethod<VariableComponentType, double>                      BossakStepMethodComponentType;
  typedef SimoStepMethod<VariableComponentType, double>                          SimoStepMethodComponentType;
  typedef EmcStepMethod<VariableComponentType, double>                            EmcStepMethodComponentType;

  typedef StaticStepRotationMethod<VariableComponentType, double>      StaticStepRotationMethodComponentType;
  typedef NewmarkStepRotationMethod<VariableComponentType, double>    NewmarkStepRotationMethodComponentType;
  typedef BossakStepRotationMethod<VariableComponentType, double>      BossakStepRotationMethodComponentType;
  typedef SimoStepRotationMethod<VariableComponentType, double>          SimoStepRotationMethodComponentType;
  typedef EmcStepRotationMethod<VariableComponentType, double>            EmcStepRotationMethodComponentType;


  
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

  //Time integraton methods for vector variables
  class_<TimeIntegrationMethodVectorType, TimeIntegrationMethodVectorType::Pointer>(m,"VectorTimeIntegration")
      .def(init<>())
      .def(init<const VariableVectorType&>())
      .def(init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      .def(init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
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
         TimeIntegrationMethodVectorType>(m,"StaticVectorIntegration")
      .def(init<>())
      .def(init<const VariableVectorType&>())
      ;

  class_<NewmarkMethodVectorType, typename NewmarkMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"NewmarkVectorIntegration")
      .def(init<>())
      .def(init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      .def(init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      ;

  class_<BossakMethodVectorType, typename BossakMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"BossakVectorIntegration")
      .def(init<>())
      .def(init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      .def(init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      ;

  class_<SimoMethodVectorType, typename SimoMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"SimoVectorIntegration")
      .def(init<>())
      .def(init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      .def(init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      ;

  class_<BackwardEulerMethodVectorType, typename BackwardEulerMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"BackwardEulerVectorIntegration")
      .def(init<>())
      .def(init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      .def(init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      ;

  class_<BdfMethodVectorType, typename BdfMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"BdfVectorIntegration")
      .def(init<>())
      .def(init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      .def(init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      ;

  class_<StaticStepMethodVectorType, typename StaticStepMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"StaticStepVectorIntegration")
      .def(init<>())
      .def(init<const VariableVectorType&>())
      ;

  class_<NewmarkStepMethodVectorType, typename NewmarkStepMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"NewmarkStepVectorIntegration")
      .def(init<>())
      .def(init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      .def(init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      ;

  class_<BossakStepMethodVectorType, typename BossakStepMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"BossakStepVectorIntegration")
      .def(init<>())
      .def(init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      .def(init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      ;

  class_<SimoStepMethodVectorType, typename SimoStepMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"SimoStepVectorIntegration")
      .def(init<>())
      .def(init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      .def(init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      ;

  class_<EmcStepMethodVectorType, typename EmcStepMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"EmcStepVectorIntegration")
      .def(init<>())
      .def(init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      .def(init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      ;

  class_<StaticStepRotationMethodVectorType, typename StaticStepRotationMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"StaticStepRotationVectorIntegration")
      .def(init<>())
      .def(init<const VariableVectorType&>())
      ;

  class_<NewmarkStepRotationMethodVectorType, typename NewmarkStepRotationMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"NewmarkStepRotationVectorIntegration")
      .def(init<>())
      .def(init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      .def(init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      ;

  class_<BossakStepRotationMethodVectorType, typename BossakStepRotationMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"BossakStepRotationVectorIntegration")
      .def(init<>())
      .def(init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      .def(init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      ;

  class_<SimoStepRotationMethodVectorType, typename SimoStepRotationMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"SimoStepRotationVectorIntegration")
      .def(init<>())
      .def(init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      .def(init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      ;

  class_<EmcStepRotationMethodVectorType, typename EmcStepRotationMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"EmcStepRotationVectorIntegration")
      .def(init<>())
      .def(init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      .def(init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      ;
  
  //Time integration methods for scalar variables
  class_<TimeIntegrationMethodScalarType, TimeIntegrationMethodScalarType::Pointer>(m,"ScalarTimeIntegration")
      .def(init<>())
      .def(init<const VariableScalarType&>())
      .def(init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      .def(init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
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
         TimeIntegrationMethodScalarType>(m,"StaticScalarIntegration")
      .def(init<>())
      .def(init<const VariableScalarType&>())
      ;

  class_<NewmarkMethodScalarType, typename NewmarkMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"NewmarkScalarIntegration")
      .def(init<>())
      .def(init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      .def(init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      ;

  class_<BossakMethodScalarType, typename BossakMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"BossakScalarIntegration")
      .def(init<>())
      .def(init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      .def(init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      ;

  class_<SimoMethodScalarType, typename SimoMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"SimoScalarIntegration")
      .def(init<>())
      .def(init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      .def(init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      ;

  class_<BackwardEulerMethodScalarType, typename BackwardEulerMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"BackwardEulerScalarIntegration")
      .def(init<>())
      .def(init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      .def(init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      ;

  class_<BdfMethodScalarType, typename BdfMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"BdfScalarIntegration")
      .def(init<>())
      .def(init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      .def(init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      ;

  class_<StaticStepMethodScalarType, typename StaticStepMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"StaticStepScalarIntegration")
      .def(init<>())
      .def(init<const VariableScalarType&>())
      ;

  class_<NewmarkStepMethodScalarType, typename NewmarkStepMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"NewmarkStepScalarIntegration")
      .def(init<>())
      .def(init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      .def(init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      ;

  class_<BossakStepMethodScalarType, typename BossakStepMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"BossakStepScalarIntegration")
      .def(init<>())
      .def(init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      .def(init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      ;

  class_<SimoStepMethodScalarType, typename SimoStepMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"SimoStepScalarIntegration")
      .def(init<>())
      .def(init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      .def(init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      ;

  class_<EmcStepMethodScalarType, typename EmcStepMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"EmcStepScalarIntegration")
      .def(init<>())
      .def(init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      .def(init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      ;

  class_<StaticStepRotationMethodScalarType, typename StaticStepRotationMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"StaticStepRotationScalarIntegration")
      //.def(init<>())
      .def(init<const VariableScalarType&>())
      ;

  // class_<NewmarkStepRotationMethodScalarType, typename NewmarkStepRotationMethodScalarType::Pointer,
  //        TimeIntegrationMethodScalarType>(m,"NewmarkStepRotationScalarIntegration")
  //     .def(init<>())
  //     .def(init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
  //     .def(init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
  //     ;

  // class_<BossakStepRotationMethodScalarType, typename BossakStepRotationMethodScalarType::Pointer,
  //        TimeIntegrationMethodScalarType>(m,"BossakStepRotationScalarIntegration")
  //     .def(init<>())
  //     .def(init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
  //     .def(init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
  //     ;

  // class_<SimoStepRotationMethodScalarType, typename SimoStepRotationMethodScalarType::Pointer,
  //        TimeIntegrationMethodScalarType>(m,"SimoStepRotationScalarIntegration")
  //     .def(init<>())
  //     .def(init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
  //     .def(init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
  //     ;

  // class_<EmcStepRotationMethodScalarType, typename EmcStepRotationMethodScalarType::Pointer,
  //        TimeIntegrationMethodScalarType>(m,"EmcStepRotationScalarIntegration")
  //     .def(init<>())     
  //     .def(init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
  //     .def(init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
  //     ;
  
  // Time integration methods for variable components
  class_<TimeIntegrationMethodComponentType, typename TimeIntegrationMethodComponentType::Pointer>(m,"ComponentTimeIntegration")
      .def(init<>())
      .def("Clone", &TimeIntegrationMethodComponentType::Clone)
      .def("SetVariable", &TimeIntegrationMethodComponentType::SetVariable)
      .def("SetFirstDerivative", &TimeIntegrationMethodComponentType::SetFirstDerivative)
      .def("SetSecondDerivative", &TimeIntegrationMethodComponentType::SetSecondDerivative)
      .def("SetVariables", &TimeIntegrationMethodComponentType::SetVariables)
      .def("SetInputVariable", &TimeIntegrationMethodComponentType::SetInputVariable)
      .def("HasStepVariable", &TimeIntegrationMethodComponentType::HasStepVariable)
      .def("SetStepVariable", &TimeIntegrationMethodComponentType::SetStepVariable)
      .def("CalculateParameters", &TimeIntegrationMethodComponentType::CalculateParameters)
      .def("SetParameters", &TimeIntegrationMethodComponentType::SetParameters)
      .def("Assign", &TimeIntegrationMethodComponentType::Assign)
      .def("Predict", &TimeIntegrationMethodComponentType::Predict)
      .def("__repr__", &TimeIntegrationMethodComponentType::Info)
      ;


  class_<StaticMethodComponentType, typename StaticMethodComponentType::Pointer,
         TimeIntegrationMethodComponentType>(m,"StaticComponentIntegration")
      .def(init<>())
      ;

  class_<NewmarkMethodComponentType, typename NewmarkMethodComponentType::Pointer,
         TimeIntegrationMethodComponentType>(m,"NewmarkComponentIntegration")
      .def(init<>())
      ;

  class_<BossakMethodComponentType, typename BossakMethodComponentType::Pointer,
         TimeIntegrationMethodComponentType>(m,"BossakComponentIntegration")
      .def(init<>())
      ;

  class_<SimoMethodComponentType, typename SimoMethodComponentType::Pointer,
         TimeIntegrationMethodComponentType>(m,"SimoComponentIntegration")
      .def(init<>())
      ;

  class_<BackwardEulerMethodComponentType, typename BackwardEulerMethodComponentType::Pointer,
         TimeIntegrationMethodComponentType>(m,"BackwardEulerComponentIntegration")
      .def(init<>())
      ;

  class_<BdfMethodComponentType, typename BdfMethodComponentType::Pointer,
         TimeIntegrationMethodComponentType>(m,"BdfComponentIntegration")
      .def(init<>())
      ;

  class_<StaticStepMethodComponentType, typename StaticStepMethodComponentType::Pointer,
         TimeIntegrationMethodComponentType>(m,"StaticStepComponentIntegration")
      .def(init<>())
      ;

  class_<NewmarkStepMethodComponentType, typename NewmarkStepMethodComponentType::Pointer,
         TimeIntegrationMethodComponentType>(m,"NewmarkStepComponentIntegration")
      .def(init<>())
      ;

  class_<BossakStepMethodComponentType, typename BossakStepMethodComponentType::Pointer,
         TimeIntegrationMethodComponentType>(m,"BossakStepComponentIntegration")
      .def(init<>())
      ;

  class_<SimoStepMethodComponentType, typename SimoStepMethodComponentType::Pointer,
         TimeIntegrationMethodComponentType>(m,"SimoStepComponentIntegration")
      .def(init<>())
      ;

  class_<EmcStepMethodComponentType, typename EmcStepMethodComponentType::Pointer,
         TimeIntegrationMethodComponentType>(m,"EmcStepComponentIntegration")
      .def(init<>())
      ;

  class_<StaticStepRotationMethodComponentType, typename StaticStepRotationMethodComponentType::Pointer,
         TimeIntegrationMethodComponentType>(m,"StaticStepRotationComponentIntegration")
      .def(init<>())
      ;

  class_<NewmarkStepRotationMethodComponentType, typename NewmarkStepRotationMethodComponentType::Pointer,
         TimeIntegrationMethodComponentType>(m,"NewmarkStepRotationComponentIntegration")
      .def(init<>())
      ;

  class_<BossakStepRotationMethodComponentType, typename BossakStepRotationMethodComponentType::Pointer,
         TimeIntegrationMethodComponentType>(m,"BossakStepRotationComponentIntegration")
      .def(init<>())
      ;

  class_<SimoStepRotationMethodComponentType, typename SimoStepRotationMethodComponentType::Pointer,
         TimeIntegrationMethodComponentType>(m,"SimoStepRotationComponentIntegration")
      .def(init<>())
      ;

  class_<EmcStepRotationMethodComponentType, typename EmcStepRotationMethodComponentType::Pointer,
         TimeIntegrationMethodComponentType>(m,"EmcStepRotationComponentIntegration")
      .def(init<>())
      ;

  // Time vector integration methods container type
  class_<VectorTimeIntegrationContainerType, VectorTimeIntegrationContainerPointerType>(m,"VectorTimeIntegrationMethods")
      .def(init<>())
      .def("Set", &VectorTimeIntegrationContainerType::Set)
      .def("Get", &VectorTimeIntegrationContainerType::Get)
      .def("Has", &VectorTimeIntegrationContainerType::Has)
      .def("__repr__", &VectorTimeIntegrationContainerType::Info)
      DECLARE_HAS_THIS_TYPE_PROCESS_INFO_PYTHON_AS_POINTER(VectorTimeIntegrationContainerType)
      DECLARE_ADD_THIS_TYPE_TO_PROCESS_INFO_PYTHON_AS_POINTER(VectorTimeIntegrationContainerType)
      DECLARE_GET_THIS_TYPE_FROM_PROCESS_INFO_PYTHON_AS_POINTER(VectorTimeIntegrationContainerType)
      ;

  //to define it as a variable
  class_<Variable<VectorTimeIntegrationContainerPointerType>, VariableData>(m,"VectorTimeIntegrationMethodsVariable")
      ;
  
  // Time vector integration methods container type
  class_<ComponentTimeIntegrationContainerType, ComponentTimeIntegrationContainerPointerType>(m,"ComponentTimeIntegrationMethods")
      .def(init<>())
      .def("Set", &ComponentTimeIntegrationContainerType::Set)
      .def("Get", &ComponentTimeIntegrationContainerType::Get)
      .def("Has", &ComponentTimeIntegrationContainerType::Has)
      .def("__repr__", &ComponentTimeIntegrationContainerType::Info)
      DECLARE_HAS_THIS_TYPE_PROCESS_INFO_PYTHON_AS_POINTER(ComponentTimeIntegrationContainerType)
      DECLARE_ADD_THIS_TYPE_TO_PROCESS_INFO_PYTHON_AS_POINTER(ComponentTimeIntegrationContainerType)
      DECLARE_GET_THIS_TYPE_FROM_PROCESS_INFO_PYTHON_AS_POINTER(ComponentTimeIntegrationContainerType)
      ;

  //to define it as a variable
  class_<Variable<ComponentTimeIntegrationContainerPointerType>, VariableData>(m,"ComponentTimeIntegrationMethodsVariable")
      ;
  
  // Time vector integration methods container type
  class_<ScalarTimeIntegrationContainerType, ScalarTimeIntegrationContainerPointerType>(m,"ScalarTimeIntegrationMethods")
      .def(init<>())
      .def("Set", &ScalarTimeIntegrationContainerType::Set)
      .def("Get", &ScalarTimeIntegrationContainerType::Get)
      .def("Has", &ScalarTimeIntegrationContainerType::Has)
      .def("__repr__", &ScalarTimeIntegrationContainerType::Info)
      DECLARE_HAS_THIS_TYPE_PROCESS_INFO_PYTHON_AS_POINTER(ScalarTimeIntegrationContainerType)
      DECLARE_ADD_THIS_TYPE_TO_PROCESS_INFO_PYTHON_AS_POINTER(ScalarTimeIntegrationContainerType)
      DECLARE_GET_THIS_TYPE_FROM_PROCESS_INFO_PYTHON_AS_POINTER(ScalarTimeIntegrationContainerType)
      ;
  
  //to define it as a variable
  class_<Variable<ScalarTimeIntegrationContainerPointerType>, VariableData>(m,"ScalarTimeIntegrationMethodsVariable")
      ;
  
}

}  // namespace Python.

} // Namespace Kratos
