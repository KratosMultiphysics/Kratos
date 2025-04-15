//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if defined(KRATOS_PYTHON)

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "custom_python/add_trilinos_strategies_to_python.h"

// Project includes
#include "trilinos_space.h"
#include "spaces/ublas_space.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/implicit_solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

// //Builder And Solver
#include "custom_strategies/builder_and_solvers/trilinos_block_builder_and_solver.h"
#include "custom_strategies/builder_and_solvers/trilinos_block_builder_and_solver_periodic.h"
#include "custom_strategies/builder_and_solvers/trilinos_elimination_builder_and_solver.h"
#include "custom_strategies/convergencecriterias/trilinos_displacement_criteria.h"

//configuration files
#include "linear_solvers/linear_solver.h"

namespace Kratos::Python
{
namespace py = pybind11;

void AddStrategies(pybind11::module& m)
{
    using TrilinosSparseSpaceType = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>;
    using TrilinosLocalSpaceType = UblasSpace<double, Matrix, Vector>;
    using TrilinosLinearSolverType = LinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType>;
    using TrilinosConvergenceCriteria = ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType>;

    using TrilinosBaseSolvingStrategyType = SolvingStrategy<TrilinosSparseSpaceType, TrilinosLocalSpaceType>;
    using TrilinosImplicitSolvingStrategyType = ImplicitSolvingStrategy<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType>;
    using TrilinosBaseSchemeType = Scheme<TrilinosSparseSpaceType, TrilinosLocalSpaceType>;
    using TrilinosBuilderAndSolverType = BuilderAndSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType>;

    //********************************************************************
    //********************************************************************
    //Builder and Solver

    // Builder and solver base class
    typedef typename ModelPart::DofsArrayType DofsArrayType;

    py::class_< TrilinosBuilderAndSolverType, typename TrilinosBuilderAndSolverType::Pointer >(m, "TrilinosResidualBasedBuilderAndSolver")
    .def(py::init<TrilinosLinearSolverType::Pointer> () )
    .def( "SetCalculateReactionsFlag", &TrilinosBuilderAndSolverType::SetCalculateReactionsFlag )
    .def( "GetCalculateReactionsFlag", &TrilinosBuilderAndSolverType::GetCalculateReactionsFlag )
    .def( "SetDofSetIsInitializedFlag", &TrilinosBuilderAndSolverType::SetDofSetIsInitializedFlag )
    .def( "GetDofSetIsInitializedFlag", &TrilinosBuilderAndSolverType::GetDofSetIsInitializedFlag )
    .def( "SetReshapeMatrixFlag", &TrilinosBuilderAndSolverType::SetReshapeMatrixFlag )
    .def( "GetReshapeMatrixFlag", &TrilinosBuilderAndSolverType::GetReshapeMatrixFlag )
    .def( "GetEquationSystemSize", &TrilinosBuilderAndSolverType::GetEquationSystemSize )
    .def( "BuildLHS", &TrilinosBuilderAndSolverType::BuildLHS )
    .def( "BuildRHS", &TrilinosBuilderAndSolverType::BuildRHS )
    .def( "Build", &TrilinosBuilderAndSolverType::Build )
    .def( "SystemSolve", &TrilinosBuilderAndSolverType::SystemSolve )
    .def( "BuildAndSolve", &TrilinosBuilderAndSolverType::BuildAndSolve )
    .def( "BuildRHSAndSolve", &TrilinosBuilderAndSolverType::BuildRHSAndSolve )
    .def( "ApplyDirichletConditions", &TrilinosBuilderAndSolverType::ApplyDirichletConditions )
    .def( "SetUpDofSet", &TrilinosBuilderAndSolverType::SetUpDofSet )
    .def( "GetDofSet",  [](TrilinosBuilderAndSolverType& self) -> DofsArrayType& {return self.GetDofSet();}, py::return_value_policy::reference_internal)
    .def( "SetUpSystem", &TrilinosBuilderAndSolverType::SetUpSystem )
    .def( "ResizeAndInitializeVectors", &TrilinosBuilderAndSolverType::ResizeAndInitializeVectors )
    .def( "InitializeSolutionStep", &TrilinosBuilderAndSolverType::InitializeSolutionStep )
    .def( "FinalizeSolutionStep", &TrilinosBuilderAndSolverType::FinalizeSolutionStep )
    .def( "CalculateReactions", &TrilinosBuilderAndSolverType::CalculateReactions )
    .def( "Clear", &TrilinosBuilderAndSolverType::Clear )
    .def( "SetEchoLevel", &TrilinosBuilderAndSolverType::SetEchoLevel )
    .def( "GetEchoLevel", &TrilinosBuilderAndSolverType::GetEchoLevel )
    ;

    using TrilinosResidualBasedEliminationBuilderAndSolverType = TrilinosResidualBasedEliminationBuilderAndSolver< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >;
    py::class_<TrilinosResidualBasedEliminationBuilderAndSolverType, typename TrilinosResidualBasedEliminationBuilderAndSolverType::Pointer, TrilinosBuilderAndSolverType>(m, "TrilinosEliminationBuilderAndSolver")
        .def(py::init<Epetra_MpiComm&, int, TrilinosLinearSolverType::Pointer > () )
    ;

    using TrilinosBlockBuilderAndSolverType = TrilinosBlockBuilderAndSolver< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >;
    py::class_<TrilinosBlockBuilderAndSolverType, typename TrilinosBlockBuilderAndSolverType::Pointer, TrilinosBuilderAndSolverType>(m, "TrilinosBlockBuilderAndSolver")
        .def(py::init<Epetra_MpiComm&, int, TrilinosLinearSolverType::Pointer > () )
        .def(py::init<Epetra_MpiComm&, TrilinosLinearSolverType::Pointer, Parameters > () )
    ;

    using TrilinosBlockBuilderAndSolverPeriodicType = TrilinosBlockBuilderAndSolverPeriodic< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >;
    py::class_<TrilinosBlockBuilderAndSolverPeriodicType, typename TrilinosBlockBuilderAndSolverPeriodicType::Pointer, TrilinosBlockBuilderAndSolverType> (m, "TrilinosBlockBuilderAndSolverPeriodic")
        .def(py::init<Epetra_MpiComm&, int, TrilinosLinearSolverType::Pointer, Kratos::Variable<int>& >() )
    ;

    // Strategy base class
    py::class_< TrilinosBaseSolvingStrategyType, typename TrilinosBaseSolvingStrategyType::Pointer >(m, "TrilinosSolvingStrategy")
        .def(py::init< ModelPart&, bool >())
        .def("Predict", &TrilinosBaseSolvingStrategyType::Predict)
        .def("Initialize", &TrilinosBaseSolvingStrategyType::Initialize)
        .def("Solve", &TrilinosBaseSolvingStrategyType::Solve)
        .def("IsConverged", &TrilinosBaseSolvingStrategyType::IsConverged)
        .def("CalculateOutputData", &TrilinosBaseSolvingStrategyType::CalculateOutputData)
        .def("SetEchoLevel", &TrilinosBaseSolvingStrategyType::SetEchoLevel)
        .def("GetEchoLevel", &TrilinosBaseSolvingStrategyType::GetEchoLevel)
        .def("SetMoveMeshFlag", &TrilinosBaseSolvingStrategyType::SetMoveMeshFlag)
        .def("MoveMeshFlag", &TrilinosBaseSolvingStrategyType::MoveMeshFlag)
        .def("MoveMesh", &TrilinosBaseSolvingStrategyType::MoveMesh)
        .def("Clear", &TrilinosBaseSolvingStrategyType::Clear)
        .def("Check", &TrilinosBaseSolvingStrategyType::Check)
        .def("InitializeSolutionStep", &TrilinosBaseSolvingStrategyType::InitializeSolutionStep)
        .def("FinalizeSolutionStep", &TrilinosBaseSolvingStrategyType::FinalizeSolutionStep)
        .def("SolveSolutionStep", &TrilinosBaseSolvingStrategyType::SolveSolutionStep)
        .def("GetModelPart", [](TrilinosBaseSolvingStrategyType& self) -> ModelPart& { return self.GetModelPart(); })
        ;

    // Implicit strategy base class
    py::class_< TrilinosImplicitSolvingStrategyType, typename TrilinosImplicitSolvingStrategyType::Pointer, TrilinosBaseSolvingStrategyType >(m, "TrilinosImplicitSolvingStrategy")
        .def("SetRebuildLevel", &TrilinosImplicitSolvingStrategyType::SetRebuildLevel)
        .def("GetRebuildLevel", &TrilinosImplicitSolvingStrategyType::GetRebuildLevel)
        .def(py::init< ModelPart&, bool >())
        ;

    using TrilinosLinearStrategy = ResidualBasedLinearStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType>;
    py::class_< TrilinosLinearStrategy , typename TrilinosLinearStrategy::Pointer, TrilinosImplicitSolvingStrategyType >
    (m,"TrilinosLinearStrategy")
        .def(py::init([](ModelPart& rModelPart, TrilinosBaseSchemeType::Pointer pScheme, TrilinosLinearSolverType::Pointer pLinearSolver, TrilinosBuilderAndSolverType::Pointer pBuilderAndSolver, bool CalculateReactionFlag, bool ReformDofSetAtEachStep, bool CalculateNormDxFlag, bool MoveMeshFlag) {
            KRATOS_WARNING("TrilinosLinearStrategy") << "Using deprecated constructor. Please use constructor without linear solver.";
            return std::shared_ptr<TrilinosLinearStrategy>(new TrilinosLinearStrategy(rModelPart, pScheme, pBuilderAndSolver, CalculateReactionFlag, ReformDofSetAtEachStep, CalculateNormDxFlag, MoveMeshFlag));
        }))
        .def(py::init< ModelPart&, TrilinosBaseSchemeType::Pointer, TrilinosBuilderAndSolverType::Pointer, bool, bool, bool, bool >())
    ;

    using TrilinosNewtonRaphsonStrategy = ResidualBasedNewtonRaphsonStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType>;
    py::class_< TrilinosNewtonRaphsonStrategy , typename TrilinosNewtonRaphsonStrategy::Pointer, TrilinosImplicitSolvingStrategyType >
    (m,"TrilinosNewtonRaphsonStrategy")
    .def(py::init([](ModelPart& rModelPart, TrilinosBaseSchemeType::Pointer pScheme, TrilinosLinearSolverType::Pointer pLinearSolver, TrilinosConvergenceCriteria::Pointer pConvergenceCriteria, TrilinosBuilderAndSolverType::Pointer pBuilderAndSolver, int MaxIterations, bool CalculateReactions, bool ReformDofSetAtEachStep, bool MoveMeshFlag) {
            KRATOS_WARNING("TrilinosNewtonRaphsonStrategy") << "Using deprecated constructor. Please use constructor without linear solver.";
            return std::shared_ptr<TrilinosNewtonRaphsonStrategy>(new TrilinosNewtonRaphsonStrategy(rModelPart, pScheme, pConvergenceCriteria, pBuilderAndSolver, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag));
        }))
    .def(py::init< ModelPart&, TrilinosBaseSchemeType::Pointer, TrilinosConvergenceCriteria::Pointer, TrilinosBuilderAndSolverType::Pointer, int, bool, bool, bool >())
    ;
}

} // namespace Kratos::Python.

#endif // KRATOS_PYTHON defined