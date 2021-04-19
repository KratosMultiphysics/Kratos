//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//

/* System includes */

/* External includes */

/* Project includes */
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_lagrange_multiplier.h"

namespace Kratos
{

///@name Type Definitions
///@{

typedef TUblasSparseSpace<double> SparseSpaceType;
typedef TUblasDenseSpace<double> LocalSpaceType;
typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;

typedef BuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> BuilderAndSolverType;
typedef ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier<SparseSpaceType,  LocalSpaceType, LinearSolverType> ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplierType;

//NOTE: here we must create persisting objects for the strategies
static ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplierType msResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier;

template<>
std::vector<Internals::RegisteredPrototypeBase<BuilderAndSolverType>> ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplierType::msPrototypes{
    Internals::RegisteredPrototype<ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplierType, BuilderAndSolverType>(ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplierType::Name(), msResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier)};

///@}

} /* namespace Kratos.*/
