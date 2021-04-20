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
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"

namespace Kratos
{

///@name Type Definitions
///@{

typedef TUblasSparseSpace<double> SparseSpaceType;
typedef TUblasDenseSpace<double> LocalSpaceType;
typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;

typedef BuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> BuilderAndSolverType;
typedef ResidualBasedEliminationBuilderAndSolver<SparseSpaceType,  LocalSpaceType, LinearSolverType> ResidualBasedEliminationBuilderAndSolverType;

//NOTE: here we must create persisting objects for the strategies
static ResidualBasedEliminationBuilderAndSolverType msResidualBasedEliminationBuilderAndSolver;

template<>
std::vector<Internals::RegisteredPrototypeBase<BuilderAndSolverType>> ResidualBasedEliminationBuilderAndSolverType::msPrototypes{
    Internals::RegisteredPrototype<ResidualBasedEliminationBuilderAndSolverType, BuilderAndSolverType>(ResidualBasedEliminationBuilderAndSolverType::Name(), msResidualBasedEliminationBuilderAndSolver)};

///@}

} /* namespace Kratos.*/
