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
#include "includes/standard_strategy_factory.h"
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

//         typedef SolvingStrategy<SpaceType,  LocalSpaceType, LinearSolverType> SolvingStrategyType;
        typedef ResidualBasedLinearStrategy<SpaceType,  LocalSpaceType, LinearSolverType> ResidualBasedLinearStrategyType;
        typedef ResidualBasedNewtonRaphsonStrategy<SpaceType,  LocalSpaceType, LinearSolverType> ResidualBasedNewtonRaphsonStrategyType;
        typedef AdaptiveResidualBasedNewtonRaphsonStrategy<SpaceType,  LocalSpaceType, LinearSolverType> AdaptiveResidualBasedNewtonRaphsonStrategyType;
        typedef LineSearchStrategy<SpaceType,  LocalSpaceType, LinearSolverType> LineSearchStrategyType;
        typedef ExplicitStrategy<SpaceType,  LocalSpaceType, LinearSolverType> ExplicitStrategyType;

        //NOTE: here we must create persisting objects for the builder and solvers
//         static auto SolvingStrategyFactory = StandardStrategyFactory<SpaceType,LocalSpaceType, LinearSolverType,SolvingStrategyType>();
        static auto ResidualBasedLinearStrategyFactory = StandardStrategyFactory<SpaceType,LocalSpaceType, LinearSolverType, ResidualBasedLinearStrategyType>();
        static auto ResidualBasedNewtonRaphsonStrategyFactory = StandardStrategyFactory<SpaceType,LocalSpaceType, LinearSolverType, ResidualBasedNewtonRaphsonStrategyType>();
        static auto AdaptiveResidualBasedNewtonRaphsonStrategyFactory = StandardStrategyFactory<SpaceType,LocalSpaceType, LinearSolverType, AdaptiveResidualBasedNewtonRaphsonStrategyType>();
        static auto LineSearchStrategyFactory = StandardStrategyFactory<SpaceType,LocalSpaceType, LinearSolverType, LineSearchStrategyType>();
        static auto ExplicitStrategyFactory = StandardStrategyFactory<SpaceType,LocalSpaceType, LinearSolverType, ExplicitStrategyType>();

        // Registration of convergence solvers
//         KRATOS_REGISTER_STRATEGY("Strategy", StrategyFactory);
        KRATOS_REGISTER_STRATEGY("ResidualBasedLinearStrategy", ResidualBasedLinearStrategyFactory);
        KRATOS_REGISTER_STRATEGY("linear_strategy", ResidualBasedLinearStrategyFactory);
        KRATOS_REGISTER_STRATEGY("ResidualBasedNewtonRaphsonStrategy", ResidualBasedNewtonRaphsonStrategyFactory);
        KRATOS_REGISTER_STRATEGY("newton_raphson_strategy", ResidualBasedNewtonRaphsonStrategyFactory);
        KRATOS_REGISTER_STRATEGY("AdaptiveResidualBasedNewtonRaphsonStrategy", AdaptiveResidualBasedNewtonRaphsonStrategyFactory);
        KRATOS_REGISTER_STRATEGY("adaptative_newton_raphson_strategy", AdaptiveResidualBasedNewtonRaphsonStrategyFactory);
        KRATOS_REGISTER_STRATEGY("LineSearchStrategy", LineSearchStrategyFactory);
        KRATOS_REGISTER_STRATEGY("line_search_strategy", LineSearchStrategyFactory);
        KRATOS_REGISTER_STRATEGY("ExplicitStrategy", ExplicitStrategyFactory);
        KRATOS_REGISTER_STRATEGY("explicit_strategy", ExplicitStrategyFactory);
    };
} // Namespace Kratos

