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

#if !defined(KRATOS_REGISTER_FACTORIES_H_INCLUDED )
#define  KRATOS_REGISTER_FACTORIES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "factories/factory.h"
#include "includes/kratos_components.h"
#include "linear_solvers/linear_solver.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/explicit_solving_strategy.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/explicit_builder.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

namespace Kratos
{
///@name Type Definitions
///@{

typedef TUblasSparseSpace<double> SparseSpaceType;
typedef TUblasDenseSpace<double> LocalSpaceType;
typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;

///@}

void RegisterStrategiesFactories();
void RegisterExplicitStrategiesFactories();
void RegisterBuilderAndSolversFactories();
void RegisterExplicitBuildersFactories();
void RegisterSchemesFactories();
void RegisterConvergenceCriteriasFactories();

typedef SolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> SolvingStrategyType;
typedef ExplicitSolvingStrategy<SparseSpaceType, LocalSpaceType> ExplicitSolvingStrategyType;
typedef BuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> BuilderAndSolverType;
typedef ExplicitBuilder<SparseSpaceType, LocalSpaceType> ExplicitBuilderType;
typedef Scheme<SparseSpaceType,LocalSpaceType> SchemeType;
typedef ConvergenceCriteria<SparseSpaceType,LocalSpaceType> ConvergenceCriteriaType;

KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) KratosComponents<SolvingStrategyType>;

#ifdef KRATOS_REGISTER_STRATEGY
#undef KRATOS_REGISTER_STRATEGY
#endif
#define KRATOS_REGISTER_STRATEGY(name, reference) \
    KratosComponents<SolvingStrategyType>::Add(name, reference);

KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) KratosComponents<ExplicitSolvingStrategyType>;

#ifdef KRATOS_REGISTER_EXPLICIT_STRATEGY
#undef KRATOS_REGISTER_EXPLICIT_STRATEGY
#endif
#define KRATOS_REGISTER_EXPLICIT_STRATEGY(name, reference) \
    KratosComponents<ExplicitSolvingStrategyType>::Add(name, reference);

KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) KratosComponents<BuilderAndSolverType>;

#ifdef KRATOS_REGISTER_BUILDER_AND_SOLVER
#undef KRATOS_REGISTER_BUILDER_AND_SOLVER
#endif
#define KRATOS_REGISTER_BUILDER_AND_SOLVER(name, reference) \
    KratosComponents<BuilderAndSolverType>::Add(name, reference);

KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) KratosComponents<ExplicitBuilderType>;

#ifdef KRATOS_REGISTER_EXPLICIT_BUILDER
#undef KRATOS_REGISTER_EXPLICIT_BUILDER
#endif
#define KRATOS_REGISTER_EXPLICIT_BUILDER(name, reference) \
    KratosComponents<ExplicitBuilderType>::Add(name, reference);

KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) KratosComponents<SchemeType>;

#ifdef KRATOS_REGISTER_SCHEME
#undef KRATOS_REGISTER_SCHEME
#endif
#define KRATOS_REGISTER_SCHEME(name, reference) \
    KratosComponents<SchemeType>::Add(name, reference);

KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) KratosComponents<ConvergenceCriteriaType>;

#ifdef KRATOS_REGISTER_CONVERGENCE_CRITERIA
#undef KRATOS_REGISTER_CONVERGENCE_CRITERIA
#endif
#define KRATOS_REGISTER_CONVERGENCE_CRITERIA(name, reference) \
    KratosComponents<ConvergenceCriteriaType>::Add(name, reference);

}  // namespace Kratos.

#endif // KRATOS_REGISTER_FACTORIES_H_INCLUDED  defined

