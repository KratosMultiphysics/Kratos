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

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"

// Utility includes
#include "includes/define.h"
#include "factories/factory.h"
#include "solving_strategies/strategies/explicit_solving_strategy.h"
#include "solving_strategies/builder_and_solvers/explicit_builder.h"
#include "spaces/ublas_space.h"

namespace Kratos
{
    namespace Testing
    {
        /// Tests

        // Spaces
        using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
        using LocalSpaceType = UblasSpace<double, Matrix, Vector>;

        /// The definition of the explicit builder
        using ExplicitBuilderType = ExplicitBuilder<SparseSpaceType,LocalSpaceType>;

        /// The definition of the factory
        using ExplicitBuilderFactoryType = Factory<ExplicitBuilderType>;

        /**
         * Checks if the ExplicitBuilder performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ExplicitBuilderFactory, KratosCoreFastSuite)
        {
            Parameters this_parameters = Parameters(R"({"name" : "explicit_builder"})");
            ExplicitBuilderType::Pointer p_explicit_builder = ExplicitBuilderFactoryType().Create(this_parameters);
            KRATOS_EXPECT_EQ(p_explicit_builder->Info(), "ExplicitBuilder");
        }

    } // namespace Testing
}  // namespace Kratos.

