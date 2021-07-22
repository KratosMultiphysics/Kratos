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
#include "solving_strategies/strategies/explicit_solving_strategy.h"

namespace Kratos
{

///@name Type Definitions
///@{

typedef TUblasSparseSpace<double> SparseSpaceType;
typedef TUblasDenseSpace<double> LocalSpaceType;

typedef ExplicitSolvingStrategy<SparseSpaceType, LocalSpaceType> ExplicitSolvingStrategyType;

//NOTE: here we must create persisting objects for the strategies
static ExplicitSolvingStrategyType msExplicitSolvingStrategy;

template<>
std::vector<Internals::RegisteredPrototypeBase<ExplicitSolvingStrategyType>> ExplicitSolvingStrategyType::msPrototypes{
    Internals::RegisteredPrototype<ExplicitSolvingStrategyType, ExplicitSolvingStrategyType>(ExplicitSolvingStrategyType::Name(), msExplicitSolvingStrategy)};

///@}

} /* namespace Kratos.*/
