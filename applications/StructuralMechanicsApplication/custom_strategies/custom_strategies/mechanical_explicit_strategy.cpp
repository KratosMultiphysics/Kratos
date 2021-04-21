//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
//kratos/license.txt
//
//  Main authors:    Klaus B Sautter (based on the work of JMCarbonel)
//
//

/* System includes */

/* External includes */

/* Project includes */
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "custom_strategies/custom_strategies/mechanical_explicit_strategy.hpp"

namespace Kratos
{

///@name Type Definitions
///@{

typedef TUblasSparseSpace<double> SparseSpaceType;
typedef TUblasDenseSpace<double> LocalSpaceType;
typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;

typedef SolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> SolvingStrategyType;
typedef MechanicalExplicitStrategy<SparseSpaceType,  LocalSpaceType, LinearSolverType> MechanicalExplicitStrategyType;

//NOTE: here we must create persisting objects for the strategies
static MechanicalExplicitStrategyType msMechanicalExplicitStrategy;

template<>
std::vector<Internals::RegisteredPrototypeBase<SolvingStrategyType>> MechanicalExplicitStrategyType::msPrototypes{
    Internals::RegisteredPrototype<MechanicalExplicitStrategyType, SolvingStrategyType>(MechanicalExplicitStrategyType::Name(), msMechanicalExplicitStrategy)};

///@}

} /* namespace Kratos.*/
