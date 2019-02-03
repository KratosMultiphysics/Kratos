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
            Parameters this_parameters = Parameters(R"({"name" : "DisplacementCriteria"})");
            ConvergenceCriteriaType::Pointer p_conv_criteria = ConvergenceCriteriaFactoryType().Create(this_parameters);
            KRATOS_CHECK_C_STRING_EQUAL((p_conv_criteria->Info()).c_str(), "DisplacementCriteria");
        }

        /**
         * Checks if the ResidualCriteria performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ResidualCriteriaFactory, KratosCoreFastSuite)
        {
            Parameters this_parameters = Parameters(R"({"name" : "ResidualCriteria"})");
            ConvergenceCriteriaType::Pointer p_conv_criteria = ConvergenceCriteriaFactoryType().Create(this_parameters);
            KRATOS_CHECK_C_STRING_EQUAL((p_conv_criteria->Info()).c_str(), "ResidualCriteria");
        }

        /**
         * Checks if the And_Criteria performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(And_CriteriaFactory, KratosCoreFastSuite)
        {
            Parameters this_parameters = Parameters(R"({"name" : "And_Criteria"})");
            ConvergenceCriteriaType::Pointer p_conv_criteria = ConvergenceCriteriaFactoryType().Create(this_parameters);
            KRATOS_CHECK_C_STRING_EQUAL((p_conv_criteria->Info()).c_str(), "And_Criteria");
        }

        /**
         * Checks if the Or_Criteria performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(Or_CriteriaFactory, KratosCoreFastSuite)
        {
            Parameters this_parameters = Parameters(R"({"name" : "Or_Criteria"})");
            ConvergenceCriteriaType::Pointer p_conv_criteria = ConvergenceCriteriaFactoryType().Create(this_parameters);
            KRATOS_CHECK_C_STRING_EQUAL((p_conv_criteria->Info()).c_str(), "Or_Criteria");
        }

    } // namespace Testing
}  // namespace Kratos.

