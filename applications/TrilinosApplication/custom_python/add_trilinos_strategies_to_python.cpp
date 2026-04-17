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
//                   Vicente Mataix Ferrandiz
//

#if defined(KRATOS_PYTHON)

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "custom_python/add_trilinos_strategies_to_python.h"

// Project includes
#include "trilinos_space.h"
#ifdef HAVE_TPETRA
#include "trilinos_space_experimental.h"
#endif
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

namespace {
    template<class TSparseSpace, class TLocalSpace, class TLinearSolverType>
    void RegisterStrategies(pybind11::module& m, const std::string& Prefix)
    {
        using ConvergenceCriteriaType = ConvergenceCriteria<TSparseSpace, TLocalSpace>;
        using BaseSchemeType = Scheme<TSparseSpace, TLocalSpace>;
        using BuilderAndSolverType = BuilderAndSolver<TSparseSpace, TLocalSpace, TLinearSolverType>;
        using BaseSolvingStrategyType = SolvingStrategy<TSparseSpace, TLocalSpace>;
        using ImplicitSolvingStrategyType = ImplicitSolvingStrategy<TSparseSpace, TLocalSpace, TLinearSolverType>;

        using DofsArrayType = typename ModelPart::DofsArrayType;

        py::class_<BuilderAndSolverType, typename BuilderAndSolverType::Pointer>(m, (Prefix + "ResidualBasedBuilderAndSolver").c_str())
            .def(py::init<typename TLinearSolverType::Pointer>())
            .def("SetCalculateReactionsFlag", &BuilderAndSolverType::SetCalculateReactionsFlag)
            .def("GetCalculateReactionsFlag", &BuilderAndSolverType::GetCalculateReactionsFlag)
            .def("SetDofSetIsInitializedFlag", &BuilderAndSolverType::SetDofSetIsInitializedFlag)
            .def("GetDofSetIsInitializedFlag", &BuilderAndSolverType::GetDofSetIsInitializedFlag)
            .def("SetReshapeMatrixFlag",       &BuilderAndSolverType::SetReshapeMatrixFlag)
            .def("GetReshapeMatrixFlag",       &BuilderAndSolverType::GetReshapeMatrixFlag)
            .def("GetEquationSystemSize",      &BuilderAndSolverType::GetEquationSystemSize)
            .def("BuildLHS",                   &BuilderAndSolverType::BuildLHS)
            .def("BuildRHS",                   &BuilderAndSolverType::BuildRHS)
            .def("Build",                      &BuilderAndSolverType::Build)
            .def("SystemSolve",                &BuilderAndSolverType::SystemSolve)
            .def("BuildAndSolve",              &BuilderAndSolverType::BuildAndSolve)
            .def("BuildRHSAndSolve",           &BuilderAndSolverType::BuildRHSAndSolve)
            .def("ApplyDirichletConditions",   &BuilderAndSolverType::ApplyDirichletConditions)
            .def("SetUpDofSet",                &BuilderAndSolverType::SetUpDofSet)
            .def("GetDofSet", [](BuilderAndSolverType& self) -> DofsArrayType& { return self.GetDofSet(); }, py::return_value_policy::reference_internal)
            .def("SetUpSystem",                &BuilderAndSolverType::SetUpSystem)
            .def("ResizeAndInitializeVectors", &BuilderAndSolverType::ResizeAndInitializeVectors)
            .def("InitializeSolutionStep",     &BuilderAndSolverType::InitializeSolutionStep)
            .def("FinalizeSolutionStep",       &BuilderAndSolverType::FinalizeSolutionStep)
            .def("CalculateReactions",         &BuilderAndSolverType::CalculateReactions)
            .def("Clear",                      &BuilderAndSolverType::Clear)
            .def("SetEchoLevel",               &BuilderAndSolverType::SetEchoLevel)
            .def("GetEchoLevel",               &BuilderAndSolverType::GetEchoLevel)
            ;

        py::class_<BaseSolvingStrategyType, typename BaseSolvingStrategyType::Pointer>(m, (Prefix + "SolvingStrategy").c_str())
            .def(py::init<ModelPart&, bool>())
            .def("Predict",                &BaseSolvingStrategyType::Predict)
            .def("Initialize",             &BaseSolvingStrategyType::Initialize)
            .def("Solve",                  &BaseSolvingStrategyType::Solve)
            .def("IsConverged",            &BaseSolvingStrategyType::IsConverged)
            .def("CalculateOutputData",    &BaseSolvingStrategyType::CalculateOutputData)
            .def("SetEchoLevel",           &BaseSolvingStrategyType::SetEchoLevel)
            .def("GetEchoLevel",           &BaseSolvingStrategyType::GetEchoLevel)
            .def("SetMoveMeshFlag",        &BaseSolvingStrategyType::SetMoveMeshFlag)
            .def("MoveMeshFlag",           &BaseSolvingStrategyType::MoveMeshFlag)
            .def("MoveMesh",               &BaseSolvingStrategyType::MoveMesh)
            .def("Clear",                  &BaseSolvingStrategyType::Clear)
            .def("Check",                  &BaseSolvingStrategyType::Check)
            .def("InitializeSolutionStep", &BaseSolvingStrategyType::InitializeSolutionStep)
            .def("FinalizeSolutionStep",   &BaseSolvingStrategyType::FinalizeSolutionStep)
            .def("SolveSolutionStep",      &BaseSolvingStrategyType::SolveSolutionStep)
            .def("GetModelPart", [](BaseSolvingStrategyType& self) -> ModelPart& { return self.GetModelPart(); })
            ;

        py::class_<ImplicitSolvingStrategyType, typename ImplicitSolvingStrategyType::Pointer, BaseSolvingStrategyType>(m, (Prefix + "ImplicitSolvingStrategy").c_str())
            .def(py::init<ModelPart&, bool>())
            .def("SetRebuildLevel", &ImplicitSolvingStrategyType::SetRebuildLevel)
            .def("GetRebuildLevel", &ImplicitSolvingStrategyType::GetRebuildLevel)
            ;

        using LinearStrategyType = ResidualBasedLinearStrategy<TSparseSpace, TLocalSpace, TLinearSolverType>;
        py::class_<LinearStrategyType, typename LinearStrategyType::Pointer, ImplicitSolvingStrategyType>(m, (Prefix + "LinearStrategy").c_str())
            .def(py::init([](ModelPart& rModelPart,
                             typename BaseSchemeType::Pointer pScheme,
                             typename TLinearSolverType::Pointer pLinearSolver,
                             typename BuilderAndSolverType::Pointer pBuilderAndSolver,
                             bool CalculateReactionFlag, bool ReformDofSetAtEachStep,
                             bool CalculateNormDxFlag, bool MoveMeshFlag) {
                KRATOS_WARNING("TrilinosLinearStrategy") << "Using deprecated constructor. Please use constructor without linear solver.";
                return std::shared_ptr<LinearStrategyType>(new LinearStrategyType(rModelPart, pScheme, pBuilderAndSolver, CalculateReactionFlag, ReformDofSetAtEachStep, CalculateNormDxFlag, MoveMeshFlag));
            }))
            .def(py::init<ModelPart&, typename BaseSchemeType::Pointer, typename BuilderAndSolverType::Pointer, bool, bool, bool, bool>())
            ;

        using NewtonRaphsonStrategyType = ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TLocalSpace, TLinearSolverType>;
        py::class_<NewtonRaphsonStrategyType, typename NewtonRaphsonStrategyType::Pointer, ImplicitSolvingStrategyType>(m, (Prefix + "NewtonRaphsonStrategy").c_str())
            .def(py::init([](ModelPart& rModelPart,
                             typename BaseSchemeType::Pointer pScheme,
                             typename TLinearSolverType::Pointer pLinearSolver,
                             typename ConvergenceCriteriaType::Pointer pConvergenceCriteria,
                             typename BuilderAndSolverType::Pointer pBuilderAndSolver,
                             int MaxIterations, bool CalculateReactions,
                             bool ReformDofSetAtEachStep, bool MoveMeshFlag) {
                KRATOS_WARNING("TrilinosNewtonRaphsonStrategy") << "Using deprecated constructor. Please use constructor without linear solver.";
                return std::shared_ptr<NewtonRaphsonStrategyType>(new NewtonRaphsonStrategyType(rModelPart, pScheme, pConvergenceCriteria, pBuilderAndSolver, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag));
            }))
            .def(py::init<ModelPart&, typename BaseSchemeType::Pointer, typename ConvergenceCriteriaType::Pointer, typename BuilderAndSolverType::Pointer, int, bool, bool, bool>())
            ;
    }
}

void AddStrategies(pybind11::module& m)
{
    using TrilinosSparseSpaceType = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>;
    using TrilinosLocalSpaceType = UblasSpace<double, Matrix, Vector>;
    using TrilinosLinearSolverType = LinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType>;
    using TrilinosBuilderAndSolverType = BuilderAndSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType>;

    RegisterStrategies<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType>(m, "Trilinos");

    // Epetra-specific concrete builder and solvers
    using TrilinosResidualBasedEliminationBuilderAndSolverType = TrilinosResidualBasedEliminationBuilderAndSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType>;
    py::class_<TrilinosResidualBasedEliminationBuilderAndSolverType, typename TrilinosResidualBasedEliminationBuilderAndSolverType::Pointer, TrilinosBuilderAndSolverType>(m, "TrilinosEliminationBuilderAndSolver")
        .def(py::init<Epetra_MpiComm&, int, TrilinosLinearSolverType::Pointer>())
        ;

    using TrilinosBlockBuilderAndSolverType = TrilinosBlockBuilderAndSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType>;
    py::class_<TrilinosBlockBuilderAndSolverType, typename TrilinosBlockBuilderAndSolverType::Pointer, TrilinosBuilderAndSolverType>(m, "TrilinosBlockBuilderAndSolver")
        .def(py::init<Epetra_MpiComm&, int, TrilinosLinearSolverType::Pointer>())
        .def(py::init<Epetra_MpiComm&, TrilinosLinearSolverType::Pointer, Parameters>())
        .def("IsConstantConstraints", &TrilinosBlockBuilderAndSolverType::IsConstantConstraints)
        .def("SetConstantConstraints", &TrilinosBlockBuilderAndSolverType::SetConstantConstraints)
        ;

    using TrilinosBlockBuilderAndSolverPeriodicType = TrilinosBlockBuilderAndSolverPeriodic<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType>;
    py::class_<TrilinosBlockBuilderAndSolverPeriodicType, typename TrilinosBlockBuilderAndSolverPeriodicType::Pointer, TrilinosBlockBuilderAndSolverType>(m, "TrilinosBlockBuilderAndSolverPeriodic")
        .def(py::init<Epetra_MpiComm&, int, TrilinosLinearSolverType::Pointer, Kratos::Variable<int>&>())
        ;

#ifdef HAVE_TPETRA
    using TrilinosExperimentalSparseSpaceType = TrilinosSpaceExperimental<Tpetra::FECrsMatrix<>, Tpetra::FEMultiVector<>>;
    using TrilinosExperimentalLinearSolverType = LinearSolver<TrilinosExperimentalSparseSpaceType, TrilinosLocalSpaceType>;

    RegisterStrategies<TrilinosExperimentalSparseSpaceType, TrilinosLocalSpaceType, TrilinosExperimentalLinearSolverType>(m, "TrilinosExperimental");
#endif
}

} // namespace Kratos::Python.

#endif // KRATOS_PYTHON defined