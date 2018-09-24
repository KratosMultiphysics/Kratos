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
#include "includes/global_variables.h"

// Utility includes
#include "includes/ublas_interface.h"

namespace Kratos
{
    namespace Testing
    {
        /// Tests

        /** Checks if it the alias operation work properly with bounded matrix
         */
        KRATOS_TEST_CASE_IN_SUITE(AMatrixAliasBounded, KratosCoreFastSuite)
        {
            constexpr double tolerance = 1e-6;

            BoundedMatrix<double, 2, 2> mat22 = ZeroMatrix(2, 2);
            mat22(0,0) = 1.0;
            mat22(1,1) = 1.0;

            BoundedMatrix<double, 2, 2> secondmat22 = ZeroMatrix(2, 2);
            secondmat22(0,1) = 1.0;
            secondmat22(1,1) = 1.0;

            const BoundedMatrix<double, 2, 2> result_mat22 = prod(mat22, secondmat22);
            mat22 = prod(mat22, secondmat22);

            KRATOS_CHECK_NEAR(norm_frobenius(result_mat22 - mat22), 0.0, tolerance);

            BoundedMatrix<double, 3, 3> mat33 = ZeroMatrix(3, 3);
            mat33(0,0) = 1.0;
            mat33(1,1) = 1.0;
            mat33(2,2) = 1.0;

            BoundedMatrix<double, 3, 3> secondmat33 = ZeroMatrix(3, 3);
            mat33(0,2) = 1.0;
            mat33(1,2) = 1.0;
            mat33(2,2) = 1.0;

            const BoundedMatrix<double, 3, 3> result_mat33 = prod(mat33, secondmat33);
            mat33 = prod(mat33, secondmat33);

            KRATOS_CHECK_NEAR(norm_frobenius(result_mat33 - mat33), 0.0, tolerance);
        }

        /** Checks if it the alias operation work properly with unbounded matrix
         */
        KRATOS_TEST_CASE_IN_SUITE(AMatrixAliasUbounded, KratosCoreFastSuite)
        {
            constexpr double tolerance = 1e-6;

            Matrix mat22 = ZeroMatrix(2, 2);
            mat22(0,0) = 1.0;
            mat22(1,1) = 1.0;

            Matrix secondmat22 = ZeroMatrix(2, 2);
            secondmat22(0,1) = 1.0;
            secondmat22(1,1) = 1.0;

            const Matrix result_mat22 = prod(mat22, secondmat22);
            mat22 = prod(mat22, secondmat22);

            KRATOS_CHECK_NEAR(norm_frobenius(result_mat22 - mat22), 0.0, tolerance);

            Matrix mat33 = ZeroMatrix(3, 3);
            mat33(0,0) = 1.0;
            mat33(1,1) = 1.0;
            mat33(2,2) = 1.0;

            Matrix secondmat33 = ZeroMatrix(3, 3);
            mat33(0,2) = 1.0;
            mat33(1,2) = 1.0;
            mat33(2,2) = 1.0;

            const Matrix result_mat33 = prod(mat33, secondmat33);
            mat33 = prod(mat33, secondmat33);

            KRATOS_CHECK_NEAR(norm_frobenius(result_mat33 - mat33), 0.0, tolerance);
        }

    } // namespace Testing
}  // namespace Kratos.

