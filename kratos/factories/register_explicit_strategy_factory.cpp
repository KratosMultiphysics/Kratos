//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "factories/register_factories.h"
#include "spaces/ublas_space.h"

// Strategies
#include "solving_strategies/strategies/explicit_solving_strategy.h"
#include "solving_strategies/strategies/explicit_solving_strategy_runge_kutta_4.h"

namespace Kratos
{
void RegisterExplicitStrategiesFactories()
{
    typedef TUblasSparseSpace<double> SpaceType;
    typedef TUblasDenseSpace<double> LocalSpaceType;

    typedef ExplicitSolvingStrategyRungeKutta4<SpaceType,  LocalSpaceType> ExplicitSolvingStrategyRungeKutta4Type;

    //NOTE: here we must create persisting objects for the strategies
    static ExplicitSolvingStrategyRungeKutta4Type msExplicitSolvingStrategyRungeKutta4;

    // Registration of explicit strategies
    KRATOS_REGISTER_EXPLICIT_STRATEGY(ExplicitSolvingStrategyRungeKutta4Type::Name(), msExplicitSolvingStrategyRungeKutta4);
};
} // Namespace Kratos
