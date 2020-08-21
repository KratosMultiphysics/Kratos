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
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"

// Utility includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "factories/base_factory.h"
#include "spaces/ublas_space.h"

namespace Kratos
{
    namespace Testing
    {
        /// Tests

        // Spaces
        using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
        using LocalSpaceType = UblasSpace<double, Matrix, Vector>;

        /// The definition of the convergence criteria
        using ConvergenceCriteriaType = ConvergenceCriteria<SparseSpaceType,LocalSpaceType>;

        /// The definition of the custom class
        using ConvergenceCriteriaFactoryType = BaseFactory<ConvergenceCriteriaType>;

        /**
         * Checks if the DisplacementCriteria performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(DisplacementCriteriaFactory, KratosCoreFastSuite)
        {
            Parameters this_parameters = Parameters(R"({"name" : "displacement_criteria"})");
            ConvergenceCriteriaType::Pointer p_conv_criteria = ConvergenceCriteriaFactoryType().Create(this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_conv_criteria->Info(), "DisplacementCriteria");
        }

        /**
         * Checks if the ResidualCriteria performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ResidualCriteriaFactory, KratosCoreFastSuite)
        {
            Parameters this_parameters = Parameters(R"({"name" : "residual_criteria"})");
            ConvergenceCriteriaType::Pointer p_conv_criteria = ConvergenceCriteriaFactoryType().Create(this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_conv_criteria->Info(), "ResidualCriteria");
        }

        /**
         * Checks if the And_Criteria performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(And_CriteriaFactory, KratosCoreFastSuite)
        {
            Parameters this_parameters = Parameters(R"({"name" : "and_criteria"})");
            ConvergenceCriteriaType::Pointer p_conv_criteria = ConvergenceCriteriaFactoryType().Create(this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_conv_criteria->Info(), "And_Criteria");
        }

        /**
         * Checks if the Or_Criteria performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(Or_CriteriaFactory, KratosCoreFastSuite)
        {
            Parameters this_parameters = Parameters(R"({"name" : "or_criteria"})");
            ConvergenceCriteriaType::Pointer p_conv_criteria = ConvergenceCriteriaFactoryType().Create(this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_conv_criteria->Info(), "Or_Criteria");
        }

        /**
         * Checks if the MixedGenericCriteria performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(MixedGenericCriteriaFactory, KratosCoreFastSuite)
        {
            Parameters this_parameters = Parameters(R"({"name" : "mixed_generic_criteria"})");
            ConvergenceCriteriaType::Pointer p_conv_criteria = ConvergenceCriteriaFactoryType().Create(this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_conv_criteria->Info(), "MixedGenericCriteria");
        }

    } // namespace Testing
}  // namespace Kratos.

