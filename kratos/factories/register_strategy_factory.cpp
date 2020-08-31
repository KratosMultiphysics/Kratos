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
#include "factories/base_factory.h"
#include "spaces/ublas_space.h"

// Strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy.h"
#include "solving_strategies/strategies/line_search_strategy.h"
//#include "solving_strategies/strategies/residualbased_arc_lenght_strategy.h"

namespace Kratos
{
void RegisterStrategies()
{
    typedef TUblasSparseSpace<double> SpaceType;
    typedef TUblasDenseSpace<double> LocalSpaceType;
    typedef LinearSolver<SpaceType, LocalSpaceType> LinearSolverType;

    typedef ResidualBasedLinearStrategy<SpaceType,  LocalSpaceType, LinearSolverType> ResidualBasedLinearStrategyType;
    typedef ResidualBasedNewtonRaphsonStrategy<SpaceType,  LocalSpaceType, LinearSolverType> ResidualBasedNewtonRaphsonStrategyType;
    typedef AdaptiveResidualBasedNewtonRaphsonStrategy<SpaceType,  LocalSpaceType, LinearSolverType> AdaptiveResidualBasedNewtonRaphsonStrategyType;
    typedef LineSearchStrategy<SpaceType,  LocalSpaceType, LinearSolverType> LineSearchStrategyType;

    //NOTE: here we must create persisting objects for the strategies
    static ResidualBasedLinearStrategyType msResidualBasedLinearStrategy;
    static ResidualBasedNewtonRaphsonStrategyType msResidualBasedNewtonRaphsonStrategy;
    static AdaptiveResidualBasedNewtonRaphsonStrategyType msAdaptiveResidualBasedNewtonRaphsonStrategy;
    static LineSearchStrategyType msLineSearchStrategy;

    // Registration of strategies
    KRATOS_REGISTER_STRATEGY(ResidualBasedLinearStrategyType::Name(), msResidualBasedLinearStrategy);
    KRATOS_REGISTER_STRATEGY(ResidualBasedNewtonRaphsonStrategyType::Name(), msResidualBasedNewtonRaphsonStrategy);
    KRATOS_REGISTER_STRATEGY(AdaptiveResidualBasedNewtonRaphsonStrategyType::Name(), msAdaptiveResidualBasedNewtonRaphsonStrategy);
    KRATOS_REGISTER_STRATEGY(LineSearchStrategyType::Name(), msLineSearchStrategy);
};
} // Namespace Kratos
