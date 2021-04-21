// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Quirin Aumann
//

/* System includes */

/* External includes */

/* Project includes */
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "custom_strategies/custom_strategies/harmonic_analysis_strategy.hpp"

namespace Kratos
{

///@name Type Definitions
///@{

typedef TUblasSparseSpace<double> SparseSpaceType;
typedef TUblasDenseSpace<double> LocalSpaceType;
typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;

typedef SolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> SolvingStrategyType;
typedef HarmonicAnalysisStrategy<SparseSpaceType,  LocalSpaceType, LinearSolverType> HarmonicAnalysisStrategyType;

//NOTE: here we must create persisting objects for the strategies
static HarmonicAnalysisStrategyType msHarmonicAnalysisStrategy;

template<>
std::vector<Internals::RegisteredPrototypeBase<SolvingStrategyType>> HarmonicAnalysisStrategyType::msPrototypes{
    Internals::RegisteredPrototype<HarmonicAnalysisStrategyType, SolvingStrategyType>(HarmonicAnalysisStrategyType::Name(), msHarmonicAnalysisStrategy)};

///@}

} /* namespace Kratos.*/