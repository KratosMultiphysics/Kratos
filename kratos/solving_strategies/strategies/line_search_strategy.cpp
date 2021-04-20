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
#include "solving_strategies/strategies/line_search_strategy.h"

namespace Kratos
{

///@name Type Definitions
///@{

typedef TUblasSparseSpace<double> SparseSpaceType;
typedef TUblasDenseSpace<double> LocalSpaceType;
typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;

typedef SolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> SolvingStrategyType;
typedef LineSearchStrategy<SparseSpaceType,  LocalSpaceType, LinearSolverType> LineSearchStrategyType;

//NOTE: here we must create persisting objects for the strategies
static LineSearchStrategyType msLineSearchStrategy;

template<>
std::vector<Internals::RegisteredPrototypeBase<SolvingStrategyType>> LineSearchStrategyType::msPrototypes{
    Internals::RegisteredPrototype<LineSearchStrategyType, SolvingStrategyType>(LineSearchStrategyType::Name(), msLineSearchStrategy)};

///@}

} /* namespace Kratos.*/
