// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Iñigo López Canalejo
//                   Klaus B. Sautter
//

/* System includes */

/* External includes */

/* Project includes */
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "custom_strategies/custom_strategies/formfinding_strategy.hpp"

namespace Kratos
{

///@name Type Definitions
///@{

typedef TUblasSparseSpace<double> SparseSpaceType;
typedef TUblasDenseSpace<double> LocalSpaceType;
typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;

typedef SolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> SolvingStrategyType;
typedef FormfindingStrategy<SparseSpaceType,  LocalSpaceType, LinearSolverType> FormfindingStrategyType;

//NOTE: here we must create persisting objects for the strategies
static FormfindingStrategyType msFormfindingStrategy;

template<>
std::vector<Internals::RegisteredPrototypeBase<SolvingStrategyType>> FormfindingStrategyType::msPrototypes{
    Internals::RegisteredPrototype<FormfindingStrategyType, SolvingStrategyType>(FormfindingStrategyType::Name(), msFormfindingStrategy)};

///@}

} /* namespace Kratos.*/