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
#include "factories/base_factory.h"
#include "spaces/ublas_space.h"

// Builder And Solver
#include "solving_strategies/builder_and_solvers/explicit_builder.h"

namespace Kratos
{
void RegisterExplicitBuildersFactories()
{
    typedef TUblasSparseSpace<double> SpaceType;
    typedef TUblasDenseSpace<double> LocalSpaceType;

    typedef ExplicitBuilder<SpaceType,  LocalSpaceType> ExplicitBuilderType;

    //NOTE: here we must create persisting objects for the builder and solvers
    static ExplicitBuilderType msExplicitBuilder;

    // Registration of builder and solvers
    KRATOS_REGISTER_EXPLICIT_BUILDER(ExplicitBuilderType::Name(), msExplicitBuilder);
};
} // Namespace Kratos

