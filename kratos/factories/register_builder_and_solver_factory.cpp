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
#include "factories/factory.h"
#include "spaces/ublas_space.h"

// Builder And Solver
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints.h"

namespace Kratos
{
void RegisterBuilderAndSolversFactories()
{
    typedef TUblasSparseSpace<double> SpaceType;
    typedef TUblasDenseSpace<double> LocalSpaceType;
    typedef LinearSolver<SpaceType, LocalSpaceType> LinearSolverType;

    typedef ResidualBasedEliminationBuilderAndSolver<SpaceType,  LocalSpaceType, LinearSolverType> ResidualBasedEliminationBuilderAndSolverType;
    typedef ResidualBasedEliminationBuilderAndSolverWithConstraints<SpaceType,  LocalSpaceType, LinearSolverType> ResidualBasedEliminationBuilderAndSolverWithConstraintsType;
    typedef ResidualBasedBlockBuilderAndSolver<SpaceType,  LocalSpaceType, LinearSolverType> ResidualBasedBlockBuilderAndSolverType;

    //NOTE: here we must create persisting objects for the builder and solvers
    static ResidualBasedEliminationBuilderAndSolverType msResidualBasedEliminationBuilderAndSolver;
    static ResidualBasedEliminationBuilderAndSolverWithConstraintsType msResidualBasedEliminationBuilderAndSolverWithConstraints;
    static ResidualBasedBlockBuilderAndSolverType msResidualBasedBlockBuilderAndSolver;

    // Registration of builder and solvers
    KRATOS_REGISTER_BUILDER_AND_SOLVER(ResidualBasedEliminationBuilderAndSolverType::Name(), msResidualBasedEliminationBuilderAndSolver);
    KRATOS_REGISTER_BUILDER_AND_SOLVER(ResidualBasedEliminationBuilderAndSolverWithConstraintsType::Name(), msResidualBasedEliminationBuilderAndSolverWithConstraints);
    KRATOS_REGISTER_BUILDER_AND_SOLVER(ResidualBasedBlockBuilderAndSolverType::Name(), msResidualBasedBlockBuilderAndSolver);
};
} // Namespace Kratos

