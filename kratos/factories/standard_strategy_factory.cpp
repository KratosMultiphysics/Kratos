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
#include "factories/standard_strategy_factory.h"
#include "spaces/ublas_space.h"

// Strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy.h"
#include "solving_strategies/strategies/line_search_strategy.h"
#include "solving_strategies/strategies/explicit_strategy.h"
//#include "solving_strategies/strategies/residualbased_arc_lenght_strategy.h"

namespace Kratos
{
    void RegisterStrategies()
    {
        typedef TUblasSparseSpace<double> SpaceType;
        typedef TUblasDenseSpace<double> LocalSpaceType;
        typedef LinearSolver<SpaceType, LocalSpaceType> LinearSolverType;

        typedef SolvingStrategy<SpaceType,  LocalSpaceType, LinearSolverType> SolvingStrategyType;
        typedef ResidualBasedLinearStrategy<SpaceType,  LocalSpaceType, LinearSolverType> ResidualBasedLinearStrategyType;
        typedef ResidualBasedNewtonRaphsonStrategy<SpaceType,  LocalSpaceType, LinearSolverType> ResidualBasedNewtonRaphsonStrategyType;
        typedef AdaptiveResidualBasedNewtonRaphsonStrategy<SpaceType,  LocalSpaceType, LinearSolverType> AdaptiveResidualBasedNewtonRaphsonStrategyType;
        typedef LineSearchStrategy<SpaceType,  LocalSpaceType, LinearSolverType> LineSearchStrategyType;
        typedef ExplicitStrategy<SpaceType,  LocalSpaceType, LinearSolverType> ExplicitStrategyType;

        //NOTE: here we must create persisting objects for the strategies
        static auto ResidualBasedLinearStrategyFactory = StandardStrategyFactory< SolvingStrategyType, ResidualBasedLinearStrategyType>();
        static auto ResidualBasedNewtonRaphsonStrategyFactory = StandardStrategyFactory< SolvingStrategyType, ResidualBasedNewtonRaphsonStrategyType>();
        static auto AdaptiveResidualBasedNewtonRaphsonStrategyFactory = StandardStrategyFactory< SolvingStrategyType, AdaptiveResidualBasedNewtonRaphsonStrategyType>();
        static auto LineSearchStrategyFactory = StandardStrategyFactory< SolvingStrategyType, LineSearchStrategyType>();
        static auto ExplicitStrategyFactory = StandardStrategyFactory< SolvingStrategyType, ExplicitStrategyType>();

        // Registration of strategies
        KRATOS_REGISTER_STRATEGY(ResidualBasedLinearStrategyType::Name(), ResidualBasedLinearStrategyFactory);
        KRATOS_REGISTER_STRATEGY(ResidualBasedNewtonRaphsonStrategyType::Name(), ResidualBasedNewtonRaphsonStrategyFactory);
        KRATOS_REGISTER_STRATEGY(AdaptiveResidualBasedNewtonRaphsonStrategyType::Name(), AdaptiveResidualBasedNewtonRaphsonStrategyFactory);
        KRATOS_REGISTER_STRATEGY(LineSearchStrategyType::Name(), LineSearchStrategyFactory);
        KRATOS_REGISTER_STRATEGY(ExplicitStrategyType::Name(), ExplicitStrategyFactory);
    };
} // Namespace Kratos
