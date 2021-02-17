//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:             BSD License
//                                       Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes

// Project includes
#include "includes/define_python.h"
#include "custom_python/add_trilinos_strategies_to_python.h"

//Trilinos includes
#include "Epetra_FEVector.h"

// Project includes
#include "trilinos_space.h"
#include "spaces/ublas_space.h"

//strategies
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

// //Builder And Solver
#include "custom_strategies/builder_and_solvers/trilinos_block_builder_and_solver.h"
#include "custom_strategies/builder_and_solvers/trilinos_block_builder_and_solver_periodic.h"
#include "custom_strategies/builder_and_solvers/trilinos_elimination_builder_and_solver.h"
#include "custom_strategies/convergencecriterias/trilinos_displacement_criteria.h"

//configuration files
#include "linear_solvers/linear_solver.h"

namespace Kratos
{
namespace Python
{
namespace py = pybind11;

void AddStrategies(pybind11::module& m)
{
    typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;
    typedef LinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosLinearSolverType;
    typedef ConvergenceCriteria< TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosConvergenceCriteria;

    typedef SolvingStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosBaseSolvingStrategyType;
    typedef Scheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosBaseSchemeType;
    typedef BuilderAndSolver< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosBuilderAndSolverType;

    //********************************************************************
    //********************************************************************
    //Builder and Solver

    // Builder and solver base class
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
    .def( "GetDofSet", &TrilinosBuilderAndSolverType::GetDofSet, py::return_value_policy::reference_internal )
    .def( "SetUpSystem", &TrilinosBuilderAndSolverType::SetUpSystem )
    .def( "ResizeAndInitializeVectors", &TrilinosBuilderAndSolverType::ResizeAndInitializeVectors )
    .def( "InitializeSolutionStep", &TrilinosBuilderAndSolverType::InitializeSolutionStep )
    .def( "FinalizeSolutionStep", &TrilinosBuilderAndSolverType::FinalizeSolutionStep )
    .def( "CalculateReactions", &TrilinosBuilderAndSolverType::CalculateReactions )
    .def( "Clear", &TrilinosBuilderAndSolverType::Clear )
    .def( "SetEchoLevel", &TrilinosBuilderAndSolverType::SetEchoLevel )
    .def( "GetEchoLevel", &TrilinosBuilderAndSolverType::GetEchoLevel )
    ;

    typedef TrilinosResidualBasedEliminationBuilderAndSolver< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosResidualBasedEliminationBuilderAndSolverType;
    py::class_<
        TrilinosResidualBasedEliminationBuilderAndSolverType,
        typename TrilinosResidualBasedEliminationBuilderAndSolverType::Pointer,
        TrilinosBuilderAndSolverType >
    (m, "TrilinosEliminationBuilderAndSolver").def(py::init<Epetra_MpiComm&, int, TrilinosLinearSolverType::Pointer > () )
    ;

    typedef TrilinosBlockBuilderAndSolver< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosBlockBuilderAndSolverType;
    py::class_<
        TrilinosBlockBuilderAndSolverType,
        typename TrilinosBlockBuilderAndSolverType::Pointer,
        TrilinosBuilderAndSolverType  >
    (m, "TrilinosBlockBuilderAndSolver").def(py::init<Epetra_MpiComm&, int, TrilinosLinearSolverType::Pointer > () )
    ;

    py::class_<
        TrilinosBlockBuilderAndSolverPeriodic< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >,
        typename TrilinosBlockBuilderAndSolverPeriodic< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::Pointer,
        TrilinosBlockBuilderAndSolverType  >
    (m, "TrilinosBlockBuilderAndSolverPeriodic").def(py::init<Epetra_MpiComm&, int, TrilinosLinearSolverType::Pointer, Kratos::Variable<int>& >() )
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
    .def("SetRebuildLevel", &TrilinosBaseSolvingStrategyType::SetRebuildLevel)
    .def("GetRebuildLevel", &TrilinosBaseSolvingStrategyType::GetRebuildLevel)
    .def("SetMoveMeshFlag", &TrilinosBaseSolvingStrategyType::SetMoveMeshFlag)
    .def("MoveMeshFlag", &TrilinosBaseSolvingStrategyType::MoveMeshFlag)
    .def("MoveMesh", &TrilinosBaseSolvingStrategyType::MoveMesh)
    .def("Clear", &TrilinosBaseSolvingStrategyType::Clear)
    .def("Check", &TrilinosBaseSolvingStrategyType::Check)
    .def("InitializeSolutionStep", &TrilinosBaseSolvingStrategyType::InitializeSolutionStep)
    .def("FinalizeSolutionStep", &TrilinosBaseSolvingStrategyType::FinalizeSolutionStep)
    .def("SolveSolutionStep", &TrilinosBaseSolvingStrategyType::SolveSolutionStep)
    .def("GetModelPart", &TrilinosBaseSolvingStrategyType::GetModelPart)
    ;

    typedef ResidualBasedLinearStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType> TrilinosLinearStrategy;
    py::class_< TrilinosLinearStrategy , typename TrilinosLinearStrategy::Pointer, TrilinosBaseSolvingStrategyType >
    (m,"TrilinosLinearStrategy")
        .def(py::init([](ModelPart& rModelPart, TrilinosBaseSchemeType::Pointer pScheme, TrilinosLinearSolverType::Pointer pLinearSolver, TrilinosBuilderAndSolverType::Pointer pBuilderAndSolver, bool CalculateReactionFlag, bool ReformDofSetAtEachStep, bool CalculateNormDxFlag, bool MoveMeshFlag) {
            KRATOS_WARNING("TrilinosLinearStrategy") << "Using deprecated constructor. Please use constructor without linear solver.";
            return std::shared_ptr<TrilinosLinearStrategy>(new TrilinosLinearStrategy(rModelPart, pScheme, pBuilderAndSolver, CalculateReactionFlag, ReformDofSetAtEachStep, CalculateNormDxFlag, MoveMeshFlag));
        }))
        .def(py::init< ModelPart&, TrilinosBaseSchemeType::Pointer, TrilinosBuilderAndSolverType::Pointer, bool, bool, bool, bool >())
    ;

    typedef ResidualBasedNewtonRaphsonStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType> TrilinosNewtonRaphsonStrategy;
    py::class_< TrilinosNewtonRaphsonStrategy , typename TrilinosNewtonRaphsonStrategy::Pointer, TrilinosBaseSolvingStrategyType >
    (m,"TrilinosNewtonRaphsonStrategy")
    .def(py::init([](ModelPart& rModelPart, TrilinosBaseSchemeType::Pointer pScheme, TrilinosLinearSolverType::Pointer pLinearSolver, TrilinosConvergenceCriteria::Pointer pConvergenceCriteria, TrilinosBuilderAndSolverType::Pointer pBuilderAndSolver, int MaxIterations, bool CalculateReactions, bool ReformDofSetAtEachStep, bool MoveMeshFlag) {
            KRATOS_WARNING("TrilinosNewtonRaphsonStrategy") << "Using deprecated constructor. Please use constructor without linear solver.";
            return std::shared_ptr<TrilinosNewtonRaphsonStrategy>(new TrilinosNewtonRaphsonStrategy(rModelPart, pScheme, pConvergenceCriteria, pBuilderAndSolver, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag));
        }))
    .def(py::init< ModelPart&, TrilinosBaseSchemeType::Pointer, TrilinosConvergenceCriteria::Pointer, TrilinosBuilderAndSolverType::Pointer, int, bool, bool, bool >())
    ;
}

} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON defined
