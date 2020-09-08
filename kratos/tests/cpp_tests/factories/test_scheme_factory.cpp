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
#include "factories/factory.h"
#include "spaces/ublas_space.h"

namespace Kratos
{
    namespace Testing
    {
        /// Tests

        // Spaces
        using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
        using LocalSpaceType = UblasSpace<double, Matrix, Vector>;

        /// The definition of the scheme
        using SchemeType = Scheme<SparseSpaceType,LocalSpaceType>;

        /// The definition of the factory
        using SchemeFactoryType = BaseFactory<SchemeType>;

        /**
         * Checks if the ResidualBasedIncrementalUpdateStaticScheme performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ResidualBasedIncrementalUpdateStaticSchemeFactory, KratosCoreFastSuite)
        {
            Parameters this_parameters = Parameters(R"({"name" : "static_scheme"})");
            SchemeType::Pointer p_scheme = SchemeFactoryType().Create(this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_scheme->Info(), "ResidualBasedIncrementalUpdateStaticScheme");
        }

        /**
         * Checks if the ResidualBasedIncrementalUpdateStaticSchemeSlip performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ResidualBasedIncrementalUpdateStaticSchemeSlipFactory, KratosCoreFastSuite)
        {
            Parameters this_parameters = Parameters(R"({"name" : "static_slip_scheme"})");
            SchemeType::Pointer p_scheme = SchemeFactoryType().Create(this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_scheme->Info(), "ResidualBasedIncrementalUpdateStaticSchemeSlip");
        }

        /**
         * Checks if the ResidualBasedBossakDisplacementScheme performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ResidualBasedBossakDisplacementSchemeFactory, KratosCoreFastSuite)
        {
            Parameters this_parameters = Parameters(R"({"name" : "bossak_scheme"})");
            SchemeType::Pointer p_scheme = SchemeFactoryType().Create(this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_scheme->Info(), "ResidualBasedBossakDisplacementScheme");
        }

        /**
         * Checks if the ResidualBasedNewmarkDisplacementScheme performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ResidualBasedNewmarkDisplacementSchemeFactory, KratosCoreFastSuite)
        {
            Parameters this_parameters = Parameters(R"({"name" : "newmark_scheme"})");
            SchemeType::Pointer p_scheme = SchemeFactoryType().Create(this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_scheme->Info(), "ResidualBasedNewmarkDisplacementScheme");
        }

        /**
         * Checks if the ResidualBasedPseudoStaticDisplacementScheme performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ResidualBasedPseudoStaticDisplacementSchemeFactory, KratosCoreFastSuite)
        {
            Parameters this_parameters = Parameters(R"({"name" : "pseudo_static_scheme", "rayleigh_beta_variable" : "PRESSURE"})");
            SchemeType::Pointer p_scheme = SchemeFactoryType().Create(this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_scheme->Info(), "ResidualBasedPseudoStaticDisplacementScheme");
        }

        /**
         * Checks if the ResidualBasedBDFDisplacementScheme performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ResidualBasedBDFDisplacementSchemeFactory, KratosCoreFastSuite)
        {
            Parameters this_parameters = Parameters(R"({"name" : "bdf_displacement_scheme"})");
            SchemeType::Pointer p_scheme = SchemeFactoryType().Create(this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_scheme->Info(), "ResidualBasedBDFDisplacementScheme");
        }

        /**
         * Checks if the ResidualBasedBDFCustomScheme performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ResidualBasedBDFCustomSchemeFactory, KratosCoreFastSuite)
        {
            Parameters this_parameters = Parameters(R"({"name" : "bdf_scheme"})");
            SchemeType::Pointer p_scheme = SchemeFactoryType().Create(this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_scheme->Info(), "ResidualBasedBDFCustomScheme");
        }

    } // namespace Testing
}  // namespace Kratos.

