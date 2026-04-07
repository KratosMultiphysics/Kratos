//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

/* System includes */

/* External includes */

/* Project includes */
#include "spaces/ublas_space.h"
#include "solving_strategies/builder_and_solvers/explicit_builder.h"

namespace Kratos
{

///@name Type Definitions
///@{

using SparseSpaceType = TUblasSparseSpace<double>;
using LocalSpaceType = TUblasDenseSpace<double>;

using ExplicitBuilderType = ExplicitBuilder<SparseSpaceType, LocalSpaceType>;

template class KRATOS_API(KRATOS_CORE) ExplicitBuilder<SparseSpaceType, LocalSpaceType>;

//NOTE: here we must create persisting objects for the builder and solvers
static ExplicitBuilderType msExplicitBuilder;

template<>
std::vector<Internals::RegisteredPrototypeBase<ExplicitBuilderType>> ExplicitBuilderType::msPrototypes{
    Internals::RegisteredPrototype<ExplicitBuilderType, ExplicitBuilderType>(ExplicitBuilderType::Name(), msExplicitBuilder)};

///@}

} /* namespace Kratos.*/
