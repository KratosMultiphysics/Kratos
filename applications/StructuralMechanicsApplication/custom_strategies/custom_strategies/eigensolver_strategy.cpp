//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main author:    Michael Andre, https://github.com/msandre
//

/* System includes */

/* External includes */

/* Project includes */
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "custom_strategies/custom_strategies/eigensolver_strategy.hpp"

namespace Kratos
{

///@name Type Definitions
///@{

typedef TUblasSparseSpace<double> SparseSpaceType;
typedef TUblasDenseSpace<double> LocalSpaceType;
typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;

typedef SolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> SolvingStrategyType;
typedef EigensolverStrategy<SparseSpaceType,  LocalSpaceType, LinearSolverType> EigensolverStrategyType;

//NOTE: here we must create persisting objects for the strategies
static EigensolverStrategyType msEigensolverStrategy;

template<>
std::vector<Internals::RegisteredPrototypeBase<SolvingStrategyType>> EigensolverStrategyType::msPrototypes{
    Internals::RegisteredPrototype<EigensolverStrategyType, SolvingStrategyType>(EigensolverStrategyType::Name(), msEigensolverStrategy)};

///@}

} /* namespace Kratos.*/