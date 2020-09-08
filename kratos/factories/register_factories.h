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
#include "factories/base_factory.h"
#include "includes/kratos_components.h"
#include "solving_strategies/strategies/explicit_solving_strategy.h"
#include "solving_strategies/builder_and_solvers/explicit_builder.h"
#include "spaces/ublas_space.h"

namespace Kratos
{
///@name Type Definitions
///@{

typedef TUblasSparseSpace<double> SparseSpaceType;
typedef TUblasDenseSpace<double> LocalSpaceType;

///@}

void RegisterExplicitStrategiesFactories();
void RegisterExplicitBuildersFactories();

typedef ExplicitSolvingStrategy<SparseSpaceType, LocalSpaceType> ExplicitSolvingStrategyType;
typedef ExplicitBuilder<SparseSpaceType, LocalSpaceType> ExplicitBuilderType;

KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) KratosComponents<ExplicitSolvingStrategyType>;

#ifdef KRATOS_REGISTER_EXPLICIT_STRATEGY
#undef KRATOS_REGISTER_EXPLICIT_STRATEGY
#endif
#define KRATOS_REGISTER_EXPLICIT_STRATEGY(name, reference) \
    KratosComponents<ExplicitSolvingStrategyType>::Add(name, reference);

KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) KratosComponents<ExplicitBuilderType>;

#ifdef KRATOS_REGISTER_EXPLICIT_BUILDER
#undef KRATOS_REGISTER_EXPLICIT_BUILDER
#endif
#define KRATOS_REGISTER_EXPLICIT_BUILDER(name, reference) \
    KratosComponents<ExplicitBuilderType>::Add(name, reference);

}  // namespace Kratos.

#endif // KRATOS_REGISTER_FACTORIES_H_INCLUDED  defined
