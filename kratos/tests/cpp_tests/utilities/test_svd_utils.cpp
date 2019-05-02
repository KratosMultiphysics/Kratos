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
#include "utilities/svd_utils.h"

// TODO: Move thid to the proper folder

namespace Kratos
{
    namespace Testing
    {
        /// Tests

        /** Checks if it calculates the SVD of a matrix 2x2
         * Checks if it calculates the SVD of a matrix 2x2
         */
        KRATOS_TEST_CASE_IN_SUITE(SVDUtilsjacobiSVD2x2Test, KratosCoreFastSuite)
        {
            constexpr double tolerance = 1e-6;

            Matrix a_matrix, u_matrix, s_matrix, v_matrix;

            a_matrix.resize(2, 2,false);
            a_matrix(0,0) = 0.57690;
            a_matrix(0,1) = 0.28760;
            a_matrix(1,0) = 0.72886;
            a_matrix(1,1) = 0.40541;

            SVDUtils<double>::JacobiSingularValueDecomposition(a_matrix, u_matrix, s_matrix, v_matrix);

            // Check decomposition is correct
            const Matrix auxmat22 = prod(u_matrix, Matrix(prod(s_matrix,v_matrix)));

            for (std::size_t i = 0; i < 2; ++i) {
                for (std::size_t j = i; j < 2; ++j) {
                    KRATOS_CHECK_NEAR(auxmat22(i,j), a_matrix(i,j), tolerance);
                }
            }

            // Check SV are correct (value and order)
            KRATOS_CHECK_NEAR(s_matrix(0,0), 1.053846, tolerance);
            KRATOS_CHECK_NEAR(s_matrix(1,1), 0.023021, tolerance);
        }

        /** Checks if it calculates the SVD of a matrix 3x3
         * Checks if it calculates the SVD of a matrix 3x3
         */
        KRATOS_TEST_CASE_IN_SUITE(SVDUtilsJacobiSVD3x3Test, KratosCoreFastSuite)
        {
            constexpr double tolerance = 1e-4;

            Matrix a_matrix, u_matrix, s_matrix, v_matrix;

            a_matrix.resize(3, 3, false);
            a_matrix(0,0) = 0.57690;
            a_matrix(0,1) = 0.28760;
            a_matrix(0,2) = 0.63942;
            a_matrix(1,0) = 0.72886;
            a_matrix(1,1) = 0.40541;
            a_matrix(1,2) = 0.13415;
            a_matrix(2,0) = 0.81972;
            a_matrix(2,1) = 0.54501;
            a_matrix(2,2) = 0.28974;

            SVDUtils<double>::JacobiSingularValueDecomposition(a_matrix, u_matrix, s_matrix, v_matrix);

            // Check decomposition is correct
            Matrix auxmat33 = prod(u_matrix, Matrix(prod(s_matrix,v_matrix)));

            for (std::size_t i = 0; i < 3; ++i) {
                for (std::size_t j = i; j < 3; ++j) {
                    KRATOS_CHECK_NEAR(auxmat33(i,j), a_matrix(i,j), tolerance);
                }
            }

            // Check SV are correct (value and order)
            KRATOS_CHECK_NEAR(s_matrix(0,0), 1.554701, tolerance);
            KRATOS_CHECK_NEAR(s_matrix(1,1), 0.412674, tolerance);
            KRATOS_CHECK_NEAR(s_matrix(2,2), 0.059198, tolerance);
        }

        /** Checks if it calculates the condition number of a matrix
         * Checks if it calculates the condition number of a matrix
         */
        KRATOS_TEST_CASE_IN_SUITE(SVDUtilsConditionNumberTest, KratosCoreFastSuite)
        {
            constexpr double tolerance = 1e-4;

            Matrix a_matrix, u_matrix, s_matrix, v_matrix;

            a_matrix.resize(3, 3, false);
            a_matrix(0,0) = 0.57690;
            a_matrix(0,1) = 0.28760;
            a_matrix(0,2) = 0.63942;
            a_matrix(1,0) = 0.72886;
            a_matrix(1,1) = 0.40541;
            a_matrix(1,2) = 0.13415;
            a_matrix(2,0) = 0.81972;
            a_matrix(2,1) = 0.54501;
            a_matrix(2,2) = 0.28974;

            const double condition_number = SVDUtils<double>::SVDConditionNumber(a_matrix); // NOTE: Considering the default tolerance

            // Check condition number is correct
            KRATOS_CHECK_NEAR(condition_number, 26.2607, tolerance);
        }

    } // namespace Testing
}  // namespace Kratos.

