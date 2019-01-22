//
//   Project Name:        KratosSolversApplication $
//   Created by:          $Author:     JMCarbonell $
//   Last modified by:    $Co-Author:              $
//   Date:                $Date:      January 2019 $
//   Revision:            $Revision:           0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"
#include "utilities/openmp_utils.h"

#include "spaces/ublas_space.h"

// Solution strategies
#include "custom_solvers/solution_strategies/newton_raphson_strategy.hpp"
#include "custom_solvers/solution_strategies/line_search_strategy.hpp"
#include "custom_solvers/solution_strategies/explicit_strategy.hpp"
#include "custom_solvers/solution_strategies/eigensolver_strategy.hpp"
#include "custom_solvers/solution_strategies/segregated_strategy.hpp"

// Solution builders and solvers
#include "linear_system/system_builders/block_builder_and_solver.hpp"
#include "linear_system/system_builders/reduction_builder_and_solver.hpp"
#include "linear_system/system_builders/explicit_builder_and_solver.hpp"

// Convergence criteria
#include "custom_solvers/convergence_criteria/residual_criterion.hpp"
#include "custom_solvers/convergence_criteria/dofs_criterion.hpp"
#include "custom_solvers/convergence_criteria/composite_criterion.hpp"

// Solution schemes
#include "custom_solvers/solution_schemes/static_scheme.hpp"
#include "custom_solvers/solution_schemes/dynamic_scheme.hpp"
#include "custom_solvers/solution_schemes/eigensolver_scheme.hpp"

#include "custom_solvers/solution_schemes/explicit_central_differences_scheme.hpp"

// Linear solvers
#include "linear_solvers/linear_solver.h"

// Time integration method
#include "custom_solvers/time_integration_methods/backward_euler_method.hpp"
#include "custom_solvers/time_integration_methods/bdf_method.hpp"
#include "custom_solvers/time_integration_methods/simo_method.hpp"

#include "custom_solvers/time_integration_methods/static_step_rotation_method.hpp"
#include "custom_solvers/time_integration_methods/newmark_step_rotation_method.hpp"
#include "custom_solvers/time_integration_methods/bossak_step_rotation_method.hpp"
#include "custom_solvers/time_integration_methods/simo_step_rotation_method.hpp"
#include "custom_solvers/time_integration_methods/emc_step_rotation_method.hpp"

// Time integration methods container
#include "custom_solvers/time_integration_methods/time_integration_methods_container.hpp"

//to update
//#include "custom_solvers/solution_strategies/explicit_hamilton_strategy.hpp"
//#include "custom_solvers/solution_schemes/explicit_hamilton_scheme.hpp"
//#include "linear_system/system_builder/explicit_hamilton_builder_and_solver.hpp"


namespace Kratos
{

namespace Python
{
namespace py = pybind11;

//base types
typedef Kratos::Vector                                                                             DenseVectorType;
typedef Kratos::Matrix                                                                             DenseMatrixType;
typedef boost::numeric::ublas::vector<double>                                                     SparseVectorType;
typedef UblasSpace<double, CompressedMatrix, SparseVectorType>                                     SparseSpaceType;
typedef UblasSpace<double, DenseMatrixType, DenseVectorType>                                        LocalSpaceType;
typedef LinearSolver<SparseSpaceType, LocalSpaceType>                                             LinearSolverType;
typedef SolutionStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>                   SolutionStrategyType;
typedef SystemBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>       SystemBuilderAndSolverType;
typedef SolutionScheme<SparseSpaceType, LocalSpaceType>                                         SolutionSchemeType;
typedef ConvergenceCriterion<SparseSpaceType, LocalSpaceType>                             ConvergenceCriterionType;

typedef SolutionStrategyType::Pointer                                                      SolutionStrategyPointer;
typedef std::vector<SolutionStrategyType::Pointer>                                     SolutionStrategiesContainer;

typedef typename ConvergenceCriterionType::Pointer                                 ConvergenceCriterionPointerType;
typedef std::vector<ConvergenceCriterionPointerType>                                  ConvergenceCriteriaContainer;


// typedef UblasSpace<double, UblasMatrix, UblasVector >                                          UblasLocalSpaceType;
// typedef LinearSolver<SparseSpaceType, UblasLocalSpaceType>                                   UblasLinearSolverType;
// typedef SolutionStrategy<SparseSpaceType, UblasLocalSpaceType, UblasLinearSolverType>    UblasSolutionStrategyType;
// typedef SolutionScheme<SparseSpaceType, UblasLocalSpaceType>                               UblasSolutionSchemeType;

void  AddCustomStrategiesToPython(pybind11::module& m)
{

  // Solution strategy types
  typedef SegregatedStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>             SegregatedStrategyType;
  typedef LinearStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>                     LinearStrategyType;
  typedef NewtonRaphsonStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>       NewtonRaphsonStrategyType;
  typedef LineSearchSolutionStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>     LineSearchStrategyType;
  typedef ExplicitSolutionStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>         ExplicitStrategyType;
  typedef EigensolverStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>           EigensolverStrategyType;

  // Solution builder_and_solver types
  typedef ReductionBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>               ReductionBuilderAndSolverType;
  typedef BlockBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>                       BlockBuilderAndSolverType;
  typedef ExplicitBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>                 ExplicitBuilderAndSolverType;

  // Solution scheme types
  typedef ExplicitCentralDifferencesScheme<SparseSpaceType, LocalSpaceType>    ExplicitCentralDifferencesSchemeType;
  typedef EigensolverScheme<SparseSpaceType, LocalSpaceType>                                  EigensolverSchemeType;

  typedef StaticScheme<SparseSpaceType, LocalSpaceType>                                            StaticSchemeType;
  typedef DynamicScheme<SparseSpaceType, LocalSpaceType>                                          DynamicSchemeType;

  // Custom convergence criterion types
  typedef ResidualCriterion<SparseSpaceType, LocalSpaceType>                                  ResidualCriterionType;
  typedef DofsCriterion<SparseSpaceType, LocalSpaceType>                                          DofsCriterionType;
  typedef CompositeCriterion<SparseSpaceType, LocalSpaceType>                                CompositeCriterionType;

  // Time integration methods for vectors
  typedef array_1d<double, 3>                                                                      VectorType;
  typedef Variable<VectorType>                                                             VariableVectorType;
  typedef TimeIntegrationMethod<VariableVectorType, VectorType>               TimeIntegrationMethodVectorType;

  typedef std::vector<TimeIntegrationMethodVectorType::Pointer>                  TimeVectorIntegrationMethods;
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

  typedef std::vector<TimeIntegrationMethodScalarType::Pointer>                  TimeScalarIntegrationMethods;
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

  // to update
  //typedef ExplicitHamiltonStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>                 ExplicitHamiltonStrategyType;
  //typedef ExplicitHamiltonScheme<SparseSpaceType, LocalSpaceType>                                       ExplicitHamiltonSchemeType;
  //typedef ExplicitHamiltonBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> ExplicitHamiltonBuilderAndSolverType;

  // Explicit Hamilton
  // py::class_<ExplicitHamiltonSchemeType, typename ExplicitHamiltonSchemeType::Pointer, SolutionSchemeType> (m,"ExplicitHamiltonScheme")
  //  .def(py::init<double, double, double, bool>())
  //   ;

  // py::class_<ExplicitHamiltonStrategyType, BaseSolvingStrategyType>(m,"ExplicitHamiltonStrategy")
  //   .def(py::init<ModelPart&, SolutionSchemeType::Pointer, LinearSolverType::Pointer, bool, bool, bool>())
  //   ;

  // py::class_<ExplicitHamiltonBuilderAndSolverType, typename ExplicitHamiltonBuilderAndSolverType::Pointer, BuilderAndSolverType>(m,"ExplicitHamiltonBuilderAndSolver")
  //   .def(py::init<LinearSolverType::Pointer> ())
  //;



  //*********************CONVERGENCE CRITERION FLAGS*********************

  // Convergence Criteria Local Flags
  py::class_<CriterionLocalFlags>(m,"CriterionLocalFlags")
      .def(py::init<>())
      .def_readonly_static("INITIALIZED", &CriterionLocalFlags::INITIALIZED)
      .def_readonly_static("INCREMENTAL", &CriterionLocalFlags::INCREMENTAL)
      .def_readonly_static("CONVERGED", &CriterionLocalFlags::CONVERGED)
      .def_readonly_static("AND", &CriterionLocalFlags::AND)
      .def_readonly_static("OR", &CriterionLocalFlags::OR)
      ;

  //***************************SOLVER FLAGS******************************

  // Solver Local Flags
  py::class_<SolverLocalFlags>(m,"SolverLocalFlags")
      .def(py::init<>())
      .def_readonly_static("INITIALIZED", &SolverLocalFlags::INITIALIZED)
      .def_readonly_static("CONVERGED", &SolverLocalFlags::CONVERGED)
      .def_readonly_static("ADAPTIVE_SOLUTION", &SolverLocalFlags::ADAPTIVE_SOLUTION)
      .def_readonly_static("MOVE_MESH", &SolverLocalFlags::MOVE_MESH)
      .def_readonly_static("UPDATE_VARIABLES", &SolverLocalFlags::UPDATE_VARIABLES)
      .def_readonly_static("REFORM_DOFS", &SolverLocalFlags::REFORM_DOFS)
      .def_readonly_static("INCREMENTAL_SOLUTION", &SolverLocalFlags::INCREMENTAL_SOLUTION)
      .def_readonly_static("COMPUTE_REACTIONS", &SolverLocalFlags::COMPUTE_REACTIONS)
      .def_readonly_static("CONSTANT_SYSTEM_MATRIX", &SolverLocalFlags::CONSTANT_SYSTEM_MATRIX)
      .def_readonly_static("RAYLEIGH_DAMPING", &SolverLocalFlags::RAYLEIGH_DAMPING)
      .def_readonly_static("IMPLEX", &SolverLocalFlags::IMPLEX)
      ;

  //***********************TIME INTEGRATION FLAGS************************

  // Convergence Criteria Local Flags
  py::class_<TimeIntegrationLocalFlags>(m,"TimeIntegrationLocalFlags")
      .def(py::init<>())
      .def_readonly_static("PREDICT_PRIMARY_VARIABLE", &TimeIntegrationLocalFlags::PREDICT_PRIMARY_VARIABLE)
      ;

  //*************************STRATEGY CLASSES***************************

  //  Base Solution Strategy
  py::class_<SolutionStrategyType, typename SolutionStrategyType::Pointer, Flags>(m,"SolutionStrategy")
      .def(py::init<ModelPart&>())
      .def(py::init<ModelPart&, Flags&>())
      .def("Initialize", &SolutionStrategyType::Initialize)
      .def("InitializeSolutionStep", &SolutionStrategyType::InitializeSolutionStep)
      .def("FinalizeSolutionStep", &SolutionStrategyType::FinalizeSolutionStep)
      .def("SolveSolutionStep", &SolutionStrategyType::SolveSolutionStep )
      .def("Solve", &SolutionStrategyType::Solve)
      .def("Check", &SolutionStrategyType::Check)
      .def("Clear", &SolutionStrategyType::Clear)
      .def("SetOptions", &SolutionStrategyType::SetOptions)
      .def("GetOptions", &SolutionStrategyType::GetOptions, py::return_value_policy::reference_internal)
      .def("SetEchoLevel", &SolutionStrategyType::SetEchoLevel)
      .def("GetEchoLevel", &SolutionStrategyType::GetEchoLevel)
      ;

  //  Segregated Strategy
  py::class_<SegregatedStrategyType, typename SegregatedStrategyType::Pointer, SolutionStrategyType>(m,"SegregatedStrategy")
      .def(py::init<ModelPart&, Flags&>())
      .def(py::init<ModelPart&, Flags&,  SolutionStrategiesContainer&>())
      .def("AddStrategy", &SegregatedStrategyType::AddStrategy)
      .def("GetStrategy", &SegregatedStrategyType::GetStrategy)
      ;

  // Linear Strategy
  py::class_<LinearStrategyType, typename LinearStrategyType::Pointer, SolutionStrategyType>(m,"LinearStrategy")
      .def(py::init<ModelPart&, SolutionSchemeType::Pointer, SystemBuilderAndSolverType::Pointer, Flags&>())
      .def(py::init<ModelPart&, SolutionSchemeType::Pointer, LinearSolverType::Pointer, Flags&>())
      ;

  // Newton Raphson Strategy
  py::class_<NewtonRaphsonStrategyType, typename NewtonRaphsonStrategyType::Pointer, LinearStrategyType>(m,"NewtonRaphsonStrategy")
      .def(py::init<ModelPart&, SolutionSchemeType::Pointer, SystemBuilderAndSolverType::Pointer, ConvergenceCriterionType::Pointer, Flags&, unsigned int>())
      .def(py::init<ModelPart&, SolutionSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriterionType::Pointer, Flags&, unsigned int>())
      .def("SetMaxIterationNumber", &NewtonRaphsonStrategyType::SetMaxIterationNumber)
      .def("GetMaxIterationNumber", &NewtonRaphsonStrategyType::GetMaxIterationNumber)
      ;

  // Newton Raphson Line Search Strategy
  py::class_<LineSearchStrategyType, typename LineSearchStrategyType::Pointer, NewtonRaphsonStrategyType>(m,"LineSearchStrategy")
      .def(py::init<ModelPart&, SolutionSchemeType::Pointer, SystemBuilderAndSolverType::Pointer, ConvergenceCriterionType::Pointer, Flags&, unsigned int>())
      .def(py::init<ModelPart&, SolutionSchemeType::Pointer, SystemBuilderAndSolverType::Pointer, ConvergenceCriterionType::Pointer, Flags&, unsigned int, unsigned int>())
      .def(py::init<ModelPart&, SolutionSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriterionType::Pointer, Flags&, unsigned int>())
      .def(py::init<ModelPart&, SolutionSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriterionType::Pointer, Flags&, unsigned int, unsigned int>())
      ;

  // Explicit Strategy
  py::class_<ExplicitStrategyType, typename ExplicitStrategyType::Pointer, SolutionStrategyType>(m,"ExplicitStrategy")
      .def(py::init<ModelPart&, SolutionSchemeType::Pointer, Flags&>())
      ;


  // Eigensolver Strategy
  py::class_<EigensolverStrategyType, typename EigensolverStrategyType::Pointer, SolutionStrategyType>(m,"EigensolverStrategy")
      .def(py::init<ModelPart&, SolutionSchemeType::Pointer, SystemBuilderAndSolverType::Pointer, Flags&, bool>())
      ;


  //*******************BUILDER AND SOLVER CLASSES***********************

  // Base Builder and Solver
  py::class_<SystemBuilderAndSolverType, typename SystemBuilderAndSolverType::Pointer, Flags>
      (m,"SystemBuilderAndSolver")
      .def(py::init<LinearSolverType::Pointer> ())
      .def(py::init<>())
      .def("BuildLHS", &SystemBuilderAndSolverType::BuildLHS)
      .def("BuildRHS", &SystemBuilderAndSolverType::BuildRHS)
      .def("Build", &SystemBuilderAndSolverType::Build)
      .def("SystemSolve", &SystemBuilderAndSolverType::SystemSolve)
      .def("BuildAndSolve", &SystemBuilderAndSolverType::BuildAndSolve)
      .def("BuildRHSAndSolve", &SystemBuilderAndSolverType::BuildRHSAndSolve)
      .def("SetUpDofSet", &SystemBuilderAndSolverType::SetUpDofSet)
      .def("GetDofSet", &SystemBuilderAndSolverType::GetDofSet, py::return_value_policy::reference_internal)
      .def("SetUpSystem", &SystemBuilderAndSolverType::SetUpSystem)
      .def("SetUpSystemMatrices", &SystemBuilderAndSolverType::SetUpSystemMatrices)
      .def("InitializeSolutionStep", &SystemBuilderAndSolverType::InitializeSolutionStep)
      .def("FinalizeSolutionStep", &SystemBuilderAndSolverType::FinalizeSolutionStep)
      .def("CalculateReactions", &SystemBuilderAndSolverType::CalculateReactions)
      .def("GetEquationSystemSize", &SystemBuilderAndSolverType::GetEquationSystemSize)
      .def("Clear", &SystemBuilderAndSolverType::Clear)
      .def("Check", &SystemBuilderAndSolverType::Check)
      .def("SetEchoLevel", &SystemBuilderAndSolverType::SetEchoLevel)
      .def("GetEchoLevel", &SystemBuilderAndSolverType::GetEchoLevel)
      ;

  py::class_<ReductionBuilderAndSolverType, typename ReductionBuilderAndSolverType::Pointer, SystemBuilderAndSolverType>(m,"ReductionBuilderAndSolver")
      .def(py::init<LinearSolverType::Pointer>())
      ;

  py::class_<BlockBuilderAndSolverType, typename BlockBuilderAndSolverType::Pointer, SystemBuilderAndSolverType>(m,"BlockBuilderAndSolver")
      .def(py::init<LinearSolverType::Pointer>())
      ;

  py::class_<ExplicitBuilderAndSolverType, typename ExplicitBuilderAndSolverType::Pointer, SystemBuilderAndSolverType>(m,"ExplicitBuilderAndSolver")
      .def(py::init<>())
      ;


  //*************************SHCHEME CLASSES****************************

  // Base Scheme
  py::class_<SolutionSchemeType, typename SolutionSchemeType::Pointer, Flags>(m,"SolutionScheme")
      .def(py::init<>())
      .def(py::init<Flags&>())
      .def(py::init<TimeVectorIntegrationMethods&, Flags&>())
      .def(py::init<TimeVectorIntegrationMethods&>())
      .def(py::init<TimeScalarIntegrationMethods&, Flags&>())
      .def(py::init<TimeScalarIntegrationMethods&>())
      .def(py::init<TimeVectorIntegrationMethods&, TimeScalarIntegrationMethods&, Flags&>())
      .def(py::init<TimeVectorIntegrationMethods&, TimeScalarIntegrationMethods&>())
      .def("Initialize", &SolutionSchemeType::Initialize)
      .def("InitializeSolutionStep", &SolutionSchemeType::InitializeSolutionStep)
      .def("FinalizeSolutionStep", &SolutionSchemeType::FinalizeSolutionStep)
      .def("SetProcess", &SolutionSchemeType::SetProcess)
      .def("SetProcessVector", &SolutionSchemeType::SetProcessVector)
      .def("Update", &SolutionSchemeType::Update)
      .def("MoveMesh", &SolutionSchemeType::MoveMesh)
      .def("Check", &SolutionSchemeType::Check)
      ;

  // Static Scheme Type
  py::class_<StaticSchemeType, typename StaticSchemeType::Pointer, SolutionSchemeType>(m,"StaticScheme")
      .def(py::init<TimeVectorIntegrationMethods&, Flags&>())
      .def(py::init<TimeVectorIntegrationMethods&>())
      .def(py::init<TimeScalarIntegrationMethods&, Flags&>())
      .def(py::init<TimeScalarIntegrationMethods&>())
      .def(py::init<TimeVectorIntegrationMethods&, TimeScalarIntegrationMethods&, Flags&>())
      .def(py::init<TimeVectorIntegrationMethods&, TimeScalarIntegrationMethods&>())
      ;

  // Dynamic Scheme Type
  py::class_<DynamicSchemeType, typename DynamicSchemeType::Pointer, SolutionSchemeType>(m,"DynamicScheme")
      .def(py::init<TimeVectorIntegrationMethods&, Flags&>())
      .def(py::init<TimeVectorIntegrationMethods&>())
      .def(py::init<TimeScalarIntegrationMethods&, Flags&>())
      .def(py::init<TimeScalarIntegrationMethods&>())
      .def(py::init<TimeVectorIntegrationMethods&, TimeScalarIntegrationMethods&, Flags&>())
      .def(py::init<TimeVectorIntegrationMethods&, TimeScalarIntegrationMethods&>())
      ;

  // Explicit scheme: Central differences
  py::class_<ExplicitCentralDifferencesSchemeType, typename ExplicitCentralDifferencesSchemeType::Pointer, SolutionSchemeType>(m,"ExplicitCentralDifferencesScheme")
      .def(py::init<Flags& ,const double, const double, const double>())
      ;

  // Eigensolver Scheme Type
  py::class_<EigensolverSchemeType, typename EigensolverSchemeType::Pointer, SolutionSchemeType>(m,"EigensolverScheme")
      .def(py::init<>())
      .def(py::init<Flags&>())
      ;



  //*******************CONVERGENCE CRITERIA CLASSES*********************

  // Convergence Criterion base type
  py::class_<ConvergenceCriterionType, typename ConvergenceCriterionType::Pointer, Flags>(m,"ConvergenceCriterion")
      .def(py::init<>())
      .def("PreCriteria", &ConvergenceCriterionType::PreCriteria)
      .def("PostCriteria", &ConvergenceCriterionType::PostCriteria)
      .def("InitializeSolutionStep", &ConvergenceCriterionType::InitializeSolutionStep)
      .def("FinalizeSolutionStep", &ConvergenceCriterionType::FinalizeSolutionStep)
      .def("Check", &ConvergenceCriterionType::Check)
      .def("SetEchoLevel", &ConvergenceCriterionType::SetEchoLevel)
      ;

  py::class_<ResidualCriterionType, typename ResidualCriterionType::Pointer, ConvergenceCriterionType>(m,"ResidualCriterion")
      .def(py::init<double, double>())
      .def(py::init<const VariableScalarType&, double, double>())
      .def(py::init<const VariableVectorType&, double, double>())
      ;

  py::class_<DofsCriterionType, typename DofsCriterionType::Pointer, ConvergenceCriterionType>(m,"DofsCriterion")
      .def(py::init<double, double>())
      .def(py::init<const VariableScalarType&, double, double>())
      .def(py::init<const VariableVectorType&, double, double>())
      ;

  py::class_<CompositeCriterionType, typename CompositeCriterionType::Pointer, ConvergenceCriterionType >(m,"CompositeCriterion")
      .def(py::init<ConvergenceCriterionPointerType, ConvergenceCriterionPointerType>())
      .def(py::init<ConvergenceCriteriaContainer&>())
      ;

  //*******************TIME INTEGRATION METHODS*************************

  //Time integraton methods for vector variables
  py::class_<TimeIntegrationMethodVectorType, TimeIntegrationMethodVectorType::Pointer, Flags>(m,"VectorTimeIntegration")
      .def(py::init<const VariableVectorType&>())
      .def(py::init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      .def(py::init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      .def("Clone", &TimeIntegrationMethodVectorType::Clone)
      .def("SetInputVariable", &TimeIntegrationMethodVectorType::SetInputVariable)
      .def("HasStepVariable", &TimeIntegrationMethodVectorType::HasStepVariable)
      .def("SetStepVariable", &TimeIntegrationMethodVectorType::SetStepVariable)
      .def("GetVariableName", &TimeIntegrationMethodVectorType::GetVariableName)
      .def("GetPrimaryVariableName", &TimeIntegrationMethodVectorType::GetPrimaryVariableName)
      .def("CalculateParameters", &TimeIntegrationMethodVectorType::CalculateParameters)
      .def("SetParameters", &TimeIntegrationMethodVectorType::SetParameters)
      .def("SetFlags", &TimeIntegrationMethodVectorType::SetFlags)
      .def("Assign", &TimeIntegrationMethodVectorType::Assign)
      .def("Predict", &TimeIntegrationMethodVectorType::Predict)
      .def("__repr__", &TimeIntegrationMethodVectorType::Info)
      ;


  py::class_<StaticMethodVectorType, typename StaticMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"StaticVectorIntegration")
      .def(py::init<const VariableVectorType&>())
      ;

  py::class_<NewmarkMethodVectorType, typename NewmarkMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"NewmarkVectorIntegration")
      .def(py::init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      .def(py::init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      ;

  py::class_<BossakMethodVectorType, typename BossakMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"BossakVectorIntegration")
      .def(py::init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      .def(py::init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      ;

  py::class_<SimoMethodVectorType, typename SimoMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"SimoVectorIntegration")
      .def(py::init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      .def(py::init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      ;

  py::class_<BackwardEulerMethodVectorType, typename BackwardEulerMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"BackwardEulerVectorIntegration")
      .def(py::init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      .def(py::init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      ;

  py::class_<BdfMethodVectorType, typename BdfMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"BdfVectorIntegration")
      .def(py::init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      .def(py::init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      ;

  py::class_<StaticStepMethodVectorType, typename StaticStepMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"StaticStepVectorIntegration")
      .def(py::init<const VariableVectorType&>())
      ;

  py::class_<NewmarkStepMethodVectorType, typename NewmarkStepMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"NewmarkStepVectorIntegration")
      .def(py::init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      .def(py::init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      ;

  py::class_<BossakStepMethodVectorType, typename BossakStepMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"BossakStepVectorIntegration")
      .def(py::init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      .def(py::init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      ;

  py::class_<SimoStepMethodVectorType, typename SimoStepMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"SimoStepVectorIntegration")
      .def(py::init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      .def(py::init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      ;

  py::class_<EmcStepMethodVectorType, typename EmcStepMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"EmcStepVectorIntegration")
      .def(py::init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      .def(py::init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      ;

  py::class_<StaticStepRotationMethodVectorType, typename StaticStepRotationMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"StaticStepRotationVectorIntegration")
      .def(py::init<const VariableVectorType&>())
      ;

  py::class_<NewmarkStepRotationMethodVectorType, typename NewmarkStepRotationMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"NewmarkStepRotationVectorIntegration")
      .def(py::init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      .def(py::init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      ;

  py::class_<BossakStepRotationMethodVectorType, typename BossakStepRotationMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"BossakStepRotationVectorIntegration")
      .def(py::init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      .def(py::init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      ;

  py::class_<SimoStepRotationMethodVectorType, typename SimoStepRotationMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"SimoStepRotationVectorIntegration")
      .def(py::init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      .def(py::init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      ;

  py::class_<EmcStepRotationMethodVectorType, typename EmcStepRotationMethodVectorType::Pointer,
         TimeIntegrationMethodVectorType>(m,"EmcStepRotationVectorIntegration")
      .def(py::init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      .def(py::init<const VariableVectorType&, const VariableVectorType&, const VariableVectorType&, const VariableVectorType&>())
      ;

  //Time integration methods for scalar variables
  py::class_<TimeIntegrationMethodScalarType, TimeIntegrationMethodScalarType::Pointer, Flags>(m,"ScalarTimeIntegration")
      .def(py::init<const VariableScalarType&>())
      .def(py::init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      .def(py::init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      .def("Clone", &TimeIntegrationMethodScalarType::Clone)
      .def("SetInputVariable", &TimeIntegrationMethodScalarType::SetInputVariable)
      .def("HasStepVariable", &TimeIntegrationMethodScalarType::HasStepVariable)
      .def("SetStepVariable", &TimeIntegrationMethodScalarType::SetStepVariable)
      .def("GetVariableName", &TimeIntegrationMethodScalarType::GetVariableName)
      .def("GetPrimaryVariableName", &TimeIntegrationMethodScalarType::GetPrimaryVariableName)
      .def("CalculateParameters", &TimeIntegrationMethodScalarType::CalculateParameters)
      .def("SetParameters", &TimeIntegrationMethodScalarType::SetParameters)
      .def("SetFlags", &TimeIntegrationMethodScalarType::SetFlags)
      .def("Assign", &TimeIntegrationMethodScalarType::Assign)
      .def("Predict", &TimeIntegrationMethodScalarType::Predict)
      .def("__repr__", &TimeIntegrationMethodScalarType::Info)
      ;


  py::class_<StaticMethodScalarType, typename StaticMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"StaticScalarIntegration")
      .def(py::init<const VariableScalarType&>())
      ;

  py::class_<NewmarkMethodScalarType, typename NewmarkMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"NewmarkScalarIntegration")
      .def(py::init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      .def(py::init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      ;

  py::class_<BossakMethodScalarType, typename BossakMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"BossakScalarIntegration")
      .def(py::init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      .def(py::init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      ;

  py::class_<SimoMethodScalarType, typename SimoMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"SimoScalarIntegration")
      .def(py::init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      .def(py::init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      ;

  py::class_<BackwardEulerMethodScalarType, typename BackwardEulerMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"BackwardEulerScalarIntegration")
      .def(py::init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      .def(py::init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      ;

  py::class_<BdfMethodScalarType, typename BdfMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"BdfScalarIntegration")
      .def(py::init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      .def(py::init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      ;

  py::class_<StaticStepMethodScalarType, typename StaticStepMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"StaticStepScalarIntegration")
      .def(py::init<const VariableScalarType&>())
      ;

  py::class_<NewmarkStepMethodScalarType, typename NewmarkStepMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"NewmarkStepScalarIntegration")
      .def(py::init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      .def(py::init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      ;

  py::class_<BossakStepMethodScalarType, typename BossakStepMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"BossakStepScalarIntegration")
      .def(py::init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      .def(py::init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      ;

  py::class_<SimoStepMethodScalarType, typename SimoStepMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"SimoStepScalarIntegration")
      .def(py::init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      .def(py::init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      ;

  py::class_<EmcStepMethodScalarType, typename EmcStepMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"EmcStepScalarIntegration")
      .def(py::init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      .def(py::init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      ;

  py::class_<StaticStepRotationMethodScalarType, typename StaticStepRotationMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"StaticStepRotationScalarIntegration")
      .def(py::init<const VariableScalarType&>())
      ;

  py::class_<NewmarkStepRotationMethodScalarType, typename NewmarkStepRotationMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"NewmarkStepRotationScalarIntegration")
      .def(py::init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      .def(py::init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      ;

  py::class_<BossakStepRotationMethodScalarType, typename BossakStepRotationMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"BossakStepRotationScalarIntegration")
      .def(py::init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      .def(py::init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      ;

  py::class_<SimoStepRotationMethodScalarType, typename SimoStepRotationMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"SimoStepRotationScalarIntegration")
      .def(py::init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      .def(py::init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      ;

  py::class_<EmcStepRotationMethodScalarType, typename EmcStepRotationMethodScalarType::Pointer,
         TimeIntegrationMethodScalarType>(m,"EmcStepRotationScalarIntegration")
      .def(py::init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      .def(py::init<const VariableScalarType&, const VariableScalarType&, const VariableScalarType&, const VariableScalarType&>())
      ;

  // Time integration methods for variable components
  py::class_<TimeIntegrationMethodComponentType, typename TimeIntegrationMethodComponentType::Pointer, Flags>(m,"ComponentTimeIntegration")
      .def(py::init<const VariableComponentType&>())
      .def(py::init<const VariableComponentType&, const VariableComponentType&, const VariableComponentType&>())
      .def(py::init<const VariableComponentType&, const VariableComponentType&, const VariableComponentType&, const VariableComponentType&>())
      .def("Clone", &TimeIntegrationMethodComponentType::Clone)
      .def("SetInputVariable", &TimeIntegrationMethodComponentType::SetInputVariable)
      .def("HasStepVariable", &TimeIntegrationMethodComponentType::HasStepVariable)
      .def("SetStepVariable", &TimeIntegrationMethodComponentType::SetStepVariable)
      .def("GetVariableName", &TimeIntegrationMethodComponentType::GetVariableName)
      .def("GetPrimaryVariableName", &TimeIntegrationMethodComponentType::GetPrimaryVariableName)
      .def("CalculateParameters", &TimeIntegrationMethodComponentType::CalculateParameters)
      .def("SetParameters", &TimeIntegrationMethodComponentType::SetParameters)
      .def("SetFlags", &TimeIntegrationMethodComponentType::SetFlags)
      .def("Assign", &TimeIntegrationMethodComponentType::Assign)
      .def("Predict", &TimeIntegrationMethodComponentType::Predict)
      .def("__repr__", &TimeIntegrationMethodComponentType::Info)
      ;

  py::class_<StaticMethodComponentType, typename StaticMethodComponentType::Pointer,
         TimeIntegrationMethodComponentType>(m,"StaticComponentIntegration")
      .def(py::init<const VariableComponentType&>())
      ;

  py::class_<NewmarkMethodComponentType, typename NewmarkMethodComponentType::Pointer,
         TimeIntegrationMethodComponentType>(m,"NewmarkComponentIntegration")
      .def(py::init<const VariableComponentType&, const VariableComponentType&, const VariableComponentType&>())
      .def(py::init<const VariableComponentType&, const VariableComponentType&, const VariableComponentType&, const VariableComponentType&>())
      ;

  py::class_<BossakMethodComponentType, typename BossakMethodComponentType::Pointer,
         TimeIntegrationMethodComponentType>(m,"BossakComponentIntegration")
      .def(py::init<const VariableComponentType&, const VariableComponentType&, const VariableComponentType&>())
      .def(py::init<const VariableComponentType&, const VariableComponentType&, const VariableComponentType&, const VariableComponentType&>())
      ;

  py::class_<SimoMethodComponentType, typename SimoMethodComponentType::Pointer,
         TimeIntegrationMethodComponentType>(m,"SimoComponentIntegration")
      .def(py::init<const VariableComponentType&, const VariableComponentType&, const VariableComponentType&>())
      .def(py::init<const VariableComponentType&, const VariableComponentType&, const VariableComponentType&, const VariableComponentType&>())
      ;

  py::class_<BackwardEulerMethodComponentType, typename BackwardEulerMethodComponentType::Pointer,
         TimeIntegrationMethodComponentType>(m,"BackwardEulerComponentIntegration")
      .def(py::init<const VariableComponentType&, const VariableComponentType&, const VariableComponentType&>())
      .def(py::init<const VariableComponentType&, const VariableComponentType&, const VariableComponentType&, const VariableComponentType&>())
      ;

  py::class_<BdfMethodComponentType, typename BdfMethodComponentType::Pointer,
         TimeIntegrationMethodComponentType>(m,"BdfComponentIntegration")
      .def(py::init<const VariableComponentType&, const VariableComponentType&, const VariableComponentType&>())
      .def(py::init<const VariableComponentType&, const VariableComponentType&, const VariableComponentType&, const VariableComponentType&>())
      ;

  py::class_<StaticStepMethodComponentType, typename StaticStepMethodComponentType::Pointer,
         TimeIntegrationMethodComponentType>(m,"StaticStepComponentIntegration")
      .def(py::init<const VariableComponentType&>())
      ;

  py::class_<NewmarkStepMethodComponentType, typename NewmarkStepMethodComponentType::Pointer,
         TimeIntegrationMethodComponentType>(m,"NewmarkStepComponentIntegration")
      .def(py::init<const VariableComponentType&, const VariableComponentType&, const VariableComponentType&>())
      .def(py::init<const VariableComponentType&, const VariableComponentType&, const VariableComponentType&, const VariableComponentType&>())
      ;

  py::class_<BossakStepMethodComponentType, typename BossakStepMethodComponentType::Pointer,
         TimeIntegrationMethodComponentType>(m,"BossakStepComponentIntegration")
      .def(py::init<const VariableComponentType&, const VariableComponentType&, const VariableComponentType&>())
      .def(py::init<const VariableComponentType&, const VariableComponentType&, const VariableComponentType&, const VariableComponentType&>())
      ;

  py::class_<SimoStepMethodComponentType, typename SimoStepMethodComponentType::Pointer,
         TimeIntegrationMethodComponentType>(m,"SimoStepComponentIntegration")
      .def(py::init<const VariableComponentType&, const VariableComponentType&, const VariableComponentType&>())
      .def(py::init<const VariableComponentType&, const VariableComponentType&, const VariableComponentType&, const VariableComponentType&>())
      ;

  py::class_<EmcStepMethodComponentType, typename EmcStepMethodComponentType::Pointer,
         TimeIntegrationMethodComponentType>(m,"EmcStepComponentIntegration")
      .def(py::init<const VariableComponentType&, const VariableComponentType&, const VariableComponentType&>())
      .def(py::init<const VariableComponentType&, const VariableComponentType&, const VariableComponentType&, const VariableComponentType&>())
      ;

  py::class_<StaticStepRotationMethodComponentType, typename StaticStepRotationMethodComponentType::Pointer,
         TimeIntegrationMethodComponentType>(m,"StaticStepRotationComponentIntegration")
      .def(py::init<const VariableComponentType&>())
      ;

  py::class_<NewmarkStepRotationMethodComponentType, typename NewmarkStepRotationMethodComponentType::Pointer,
         TimeIntegrationMethodComponentType>(m,"NewmarkStepRotationComponentIntegration")
      .def(py::init<const VariableComponentType&, const VariableComponentType&, const VariableComponentType&>())
      .def(py::init<const VariableComponentType&, const VariableComponentType&, const VariableComponentType&, const VariableComponentType&>())
      ;

  py::class_<BossakStepRotationMethodComponentType, typename BossakStepRotationMethodComponentType::Pointer,
         TimeIntegrationMethodComponentType>(m,"BossakStepRotationComponentIntegration")
      .def(py::init<const VariableComponentType&, const VariableComponentType&, const VariableComponentType&>())
      .def(py::init<const VariableComponentType&, const VariableComponentType&, const VariableComponentType&, const VariableComponentType&>())
      ;

  py::class_<SimoStepRotationMethodComponentType, typename SimoStepRotationMethodComponentType::Pointer,
         TimeIntegrationMethodComponentType>(m,"SimoStepRotationComponentIntegration")
      .def(py::init<const VariableComponentType&, const VariableComponentType&, const VariableComponentType&>())
      .def(py::init<const VariableComponentType&, const VariableComponentType&, const VariableComponentType&, const VariableComponentType&>())
      ;

  py::class_<EmcStepRotationMethodComponentType, typename EmcStepRotationMethodComponentType::Pointer,
         TimeIntegrationMethodComponentType>(m,"EmcStepRotationComponentIntegration")
      .def(py::init<const VariableComponentType&, const VariableComponentType&, const VariableComponentType&>())
      .def(py::init<const VariableComponentType&, const VariableComponentType&, const VariableComponentType&, const VariableComponentType&>())
      ;

  // Time vector integration methods container type
  py::class_<VectorTimeIntegrationContainerType, VectorTimeIntegrationContainerPointerType>(m,"VectorTimeIntegrationMethods")
      .def(py::init<>())
      .def("Set", &VectorTimeIntegrationContainerType::Set)
      .def("Get", &VectorTimeIntegrationContainerType::Get)
      .def("Has", &VectorTimeIntegrationContainerType::Has)
      .def("GetMethodVariableName", &VectorTimeIntegrationContainerType::GetMethodVariableName)
      .def("__repr__", &VectorTimeIntegrationContainerType::Info)
      DECLARE_HAS_THIS_TYPE_PROCESS_INFO_PYTHON_AS_POINTER(VectorTimeIntegrationContainerType)
      DECLARE_ADD_THIS_TYPE_TO_PROCESS_INFO_PYTHON_AS_POINTER(VectorTimeIntegrationContainerType)
      DECLARE_GET_THIS_TYPE_FROM_PROCESS_INFO_PYTHON_AS_POINTER(VectorTimeIntegrationContainerType)
      ;

  //to define it as a variable
  py::class_<Variable<VectorTimeIntegrationContainerPointerType>, VariableData>(m,"VectorTimeIntegrationMethodsVariable")
      .def( "__repr__", &Variable<VectorTimeIntegrationContainerPointerType>::Info )
      ;

  // Time vector integration methods container type
  py::class_<ComponentTimeIntegrationContainerType, ComponentTimeIntegrationContainerPointerType>(m,"ComponentTimeIntegrationMethods")
      .def(py::init<>())
      .def("Set", &ComponentTimeIntegrationContainerType::Set)
      .def("Get", &ComponentTimeIntegrationContainerType::Get)
      .def("Has", &ComponentTimeIntegrationContainerType::Has)
      .def("GetMethodVariableName", &ComponentTimeIntegrationContainerType::GetMethodVariableName)
      .def("__repr__", &ComponentTimeIntegrationContainerType::Info)
      DECLARE_HAS_THIS_TYPE_PROCESS_INFO_PYTHON_AS_POINTER(ComponentTimeIntegrationContainerType)
      DECLARE_ADD_THIS_TYPE_TO_PROCESS_INFO_PYTHON_AS_POINTER(ComponentTimeIntegrationContainerType)
      DECLARE_GET_THIS_TYPE_FROM_PROCESS_INFO_PYTHON_AS_POINTER(ComponentTimeIntegrationContainerType)
      ;

  //to define it as a variable
  py::class_<Variable<ComponentTimeIntegrationContainerPointerType>, VariableData>(m,"ComponentTimeIntegrationMethodsVariable")
      .def( "__repr__", &Variable<ComponentTimeIntegrationContainerPointerType>::Info )
      ;

  // Time vector integration methods container type
  py::class_<ScalarTimeIntegrationContainerType, ScalarTimeIntegrationContainerPointerType>(m,"ScalarTimeIntegrationMethods")
      .def(py::init<>())
      .def("Set", &ScalarTimeIntegrationContainerType::Set)
      .def("Get", &ScalarTimeIntegrationContainerType::Get)
      .def("Has", &ScalarTimeIntegrationContainerType::Has)
      .def("GetMethodVariableName", &ScalarTimeIntegrationContainerType::GetMethodVariableName)
      .def("__repr__", &ScalarTimeIntegrationContainerType::Info)
      DECLARE_HAS_THIS_TYPE_PROCESS_INFO_PYTHON_AS_POINTER(ScalarTimeIntegrationContainerType)
      DECLARE_ADD_THIS_TYPE_TO_PROCESS_INFO_PYTHON_AS_POINTER(ScalarTimeIntegrationContainerType)
      DECLARE_GET_THIS_TYPE_FROM_PROCESS_INFO_PYTHON_AS_POINTER(ScalarTimeIntegrationContainerType)
      ;

  //to define it as a variable
  py::class_<Variable<ScalarTimeIntegrationContainerPointerType>, VariableData>(m,"ScalarTimeIntegrationMethodsVariable")
      .def( "__repr__", &Variable<ScalarTimeIntegrationContainerPointerType>::Info )
      ;

}

}  // namespace Python.

} // Namespace Kratos
