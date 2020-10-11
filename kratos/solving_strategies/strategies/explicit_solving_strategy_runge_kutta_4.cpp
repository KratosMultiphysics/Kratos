//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

/* System includes */

/* External includes */

/* Project includes */
#include "spaces/ublas_space.h"
#include "solving_strategies/strategies/explicit_solving_strategy_runge_kutta_4.h"

namespace Kratos
{

///@name Type Definitions
///@{

typedef TUblasSparseSpace<double> SparseSpaceType;
typedef TUblasDenseSpace<double> LocalSpaceType;

typedef ExplicitSolvingStrategy<SparseSpaceType, LocalSpaceType> ExplicitSolvingStrategyType;
typedef ExplicitSolvingStrategyRungeKutta4<SparseSpaceType, LocalSpaceType> ExplicitSolvingStrategyRungeKutta4Type;

//NOTE: here we must create persisting objects for the strategies
static ExplicitSolvingStrategyRungeKutta4Type msExplicitSolvingStrategyRungeKutta4;

template<>
std::vector<Internals::RegisteredPrototypeBase<ExplicitSolvingStrategyType>> ExplicitSolvingStrategyRungeKutta4Type::msPrototypes{
    Internals::RegisteredPrototype<ExplicitSolvingStrategyRungeKutta4Type, ExplicitSolvingStrategyType>(ExplicitSolvingStrategyRungeKutta4Type::Name(), msExplicitSolvingStrategyRungeKutta4)};

///@}

} /* namespace Kratos.*/
