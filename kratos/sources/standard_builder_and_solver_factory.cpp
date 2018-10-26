//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/standard_builder_and_solver_factory.h"
#include "spaces/ublas_space.h"

// Builder And Solver
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_constraints.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"

namespace Kratos
{
    void RegisterBuilderAndSolvers()
    {
        typedef TUblasSparseSpace<double> SpaceType;
        typedef TUblasDenseSpace<double> LocalSpaceType;
        typedef LinearSolver<SpaceType, LocalSpaceType> LinearSolverType;

//         typedef BuilderAndSolver<SpaceType,  LocalSpaceType, LinearSolverType> BuilderAndSolverType;
        typedef ResidualBasedEliminationBuilderAndSolver<SpaceType,  LocalSpaceType, LinearSolverType> ResidualBasedEliminationBuilderAndSolverType;
        typedef ResidualBasedBlockBuilderAndSolver<SpaceType,  LocalSpaceType, LinearSolverType> ResidualBasedBlockBuilderAndSolverType;
        typedef ResidualBasedBlockBuilderAndSolverWithConstraints<SpaceType,  LocalSpaceType, LinearSolverType> ResidualBasedBlockBuilderAndSolverWithConstraintsType;

        //NOTE: here we must create persisting objects for the builder and solvers
//         static auto BuilderAndSolverFactory = StandardBuilderAndSolverFactory<SpaceType,LocalSpaceType, LinearSolverType,BuilderAndSolverType>();
        static auto ResidualBasedEliminationBuilderAndSolverFactory = StandardBuilderAndSolverFactory<SpaceType,LocalSpaceType, LinearSolverType, ResidualBasedEliminationBuilderAndSolverType>();
        static auto ResidualBasedBlockBuilderAndSolverFactory = StandardBuilderAndSolverFactory<SpaceType,LocalSpaceType, LinearSolverType, ResidualBasedBlockBuilderAndSolverType>();
        static auto ResidualBasedBlockBuilderAndSolverWithConstraintsFactory = StandardBuilderAndSolverFactory<SpaceType,LocalSpaceType, LinearSolverType, ResidualBasedBlockBuilderAndSolverWithConstraintsType>();

        // Registration of convergence solvers
//         KRATOS_REGISTER_BUILDER_AND_SOLVER("BuilderAndSolver", BuilderAndSolverFactory);
        KRATOS_REGISTER_BUILDER_AND_SOLVER("ResidualBasedEliminationBuilderAndSolver", ResidualBasedEliminationBuilderAndSolverFactory);
        KRATOS_REGISTER_BUILDER_AND_SOLVER("elimination_builder_and_solver", ResidualBasedEliminationBuilderAndSolverFactory);
        KRATOS_REGISTER_BUILDER_AND_SOLVER("ResidualBasedBlockBuilderAndSolver", ResidualBasedBlockBuilderAndSolverFactory);
        KRATOS_REGISTER_BUILDER_AND_SOLVER("block_builder_and_solver", ResidualBasedBlockBuilderAndSolverFactory);
        KRATOS_REGISTER_BUILDER_AND_SOLVER("ResidualBasedBlockBuilderAndSolverWithConstraints", ResidualBasedBlockBuilderAndSolverWithConstraintsFactory);
        KRATOS_REGISTER_BUILDER_AND_SOLVER("block_builder_and_solver_with_constraints", ResidualBasedBlockBuilderAndSolverWithConstraintsFactory);
    };
} // Namespace Kratos

