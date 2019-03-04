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
#include "utilities/math_utils.h"

namespace Kratos
{
    namespace Testing
    {
        /// Tests

        /** Checks if the area of the triangle is calculated correctly using Heron equation.
         * Checks if the area of the triangle is calculated correctly using Heron equation.
         */

        KRATOS_TEST_CASE_IN_SUITE(MathUtilsHeron, KratosCoreFastSuite)
        {
            constexpr double tolerance = 1e-6;

            const double area = MathUtils<double>::Heron<false>(std::sqrt(2.0), 1.0, 1.0);

            KRATOS_CHECK_NEAR(area, 0.5, tolerance);
        }

        /** Checks if it gives you the absolute value of a given value
         * Checks if It gives you the absolute value of a given value
         */

        KRATOS_TEST_CASE_IN_SUITE(MathUtilsAbs, KratosCoreFastSuite)
        {
            const double absolute = MathUtils<double>::Abs(-1.0);

            KRATOS_CHECK_EQUAL(absolute, 1.0);
        }

        /** Checks if it gives you the minimum value of a given value
         * Checks if It gives you the minimum value of a given value
         */

        KRATOS_TEST_CASE_IN_SUITE(MathUtilsMin, KratosCoreFastSuite)
        {
            const double min = MathUtils<double>::Min(0.0,1.0);

            KRATOS_CHECK_EQUAL(min, 0.0);
        }

        /** Checks if it gives you the maximum value of a given value
         * Checks if It gives you the maximum value of a given value
         */

        KRATOS_TEST_CASE_IN_SUITE(MathUtilsMax, KratosCoreFastSuite)
        {
            const double max = MathUtils<double>::Max(0.0,1.0);

            KRATOS_CHECK_EQUAL(max, 1.0);
        }

        /** Checks if it calculates the determinant of a 1x1, 2x2, 3x3 and 4x4 matrix
         * Checks if it calculates the determinant of a 1x1, 2x2, 3x3 and 4x4 matrix
         */

        KRATOS_TEST_CASE_IN_SUITE(MathUtilsDetMat, KratosCoreFastSuite)
        {
            constexpr double tolerance = 1e-6;

            BoundedMatrix<double, 1, 1> mat11 = ZeroMatrix(1, 1);
            mat11(0,0) = 1.0;

            double det = MathUtils<double>::DetMat(mat11);

            KRATOS_CHECK_NEAR(det, 1.0, tolerance);

            BoundedMatrix<double, 2, 2> mat22 = ZeroMatrix(2, 2);
            mat22(0,0) = 1.0;
            mat22(1,1) = 1.0;

            det = MathUtils<double>::DetMat(mat22);

            KRATOS_CHECK_NEAR(det, 1.0, tolerance);

            BoundedMatrix<double, 3, 3> mat33 = ZeroMatrix(3, 3);
            mat33(0,0) = 1.0;
            mat33(1,1) = 1.0;
            mat33(2,2) = 1.0;

            det = MathUtils<double>::DetMat(mat33);

            KRATOS_CHECK_NEAR(det, 1.0, tolerance);

            BoundedMatrix<double, 4, 4> mat44 = ZeroMatrix(4, 4);
            mat44(0,0) = 1.0;
            mat44(1,1) = 1.0;
            mat44(2,2) = 1.0;
            mat44(3,3) = 1.0;

            det = MathUtils<double>::DetMat(mat44);

            KRATOS_CHECK_NEAR(det, 1.0, tolerance);
        }

        KRATOS_TEST_CASE_IN_SUITE(MathUtilsCofactor, KratosCoreFastSuite)
        {
            constexpr double tolerance = 1e-6;

            BoundedMatrix<double, 1, 1> mat11 = ZeroMatrix(1, 1);
            mat11(0,0) = 2.0;

            double cofactor = MathUtils<double>::Cofactor(mat11, 0, 0);

            KRATOS_CHECK_EQUAL(cofactor, 1.0);

            BoundedMatrix<double, 2, 2> mat22 = ZeroMatrix(2, 2);
            mat22(0,0) = -2.0; mat22(0,1) = 2.0;
            mat22(1,0) = -1.0; mat22(1,1) = 1.0;

            cofactor = MathUtils<double>::Cofactor(mat22, 1, 1);
            KRATOS_CHECK_EQUAL(cofactor, -2.0);

            cofactor = MathUtils<double>::Cofactor(mat22, 0, 1);
            KRATOS_CHECK_EQUAL(cofactor, 1.0);

            BoundedMatrix<double, 3, 3> mat33 = ZeroMatrix(3, 3);
            mat33(0,0) = -2.0; mat33(0,1) = 2.0; mat33(0,2) = -3.0;
            mat33(1,0) = -1.0; mat33(1,1) = 1.0; mat33(1,2) = 3.0;
            mat33(2,0) = 2.0; mat33(2,1) = 0.0; mat33(2,2) = -1.0;

            cofactor = MathUtils<double>::Cofactor(mat33, 2, 1);
            KRATOS_CHECK_NEAR(cofactor, 9.0, tolerance);
        }

        KRATOS_TEST_CASE_IN_SUITE(MathUtilsCofactorMatrix, KratosCoreFastSuite)
        {
            constexpr double tolerance = 1e-6;

            BoundedMatrix<double, 3, 3> mat33 = ZeroMatrix(3, 3);
            mat33(0,0) = 2.0; mat33(0,1) = 0.0; mat33(0,2) = 2.0;
            mat33(1,0) = 2.0; mat33(1,1) = 0.0; mat33(1,2) =-2.0;
            mat33(2,0) = 0.0; mat33(2,1) = 1.0; mat33(2,2) = 1.0;

            BoundedMatrix<double, 3, 3> ref33 = ZeroMatrix(3, 3);
            ref33(0,0) = 2.0; ref33(0,1) =-2.0; ref33(0,2) = 2.0;
            ref33(1,0) = 2.0; ref33(1,1) = 2.0; ref33(1,2) =-2.0;
            ref33(2,0) = 0.0; ref33(2,1) = 8.0; ref33(2,2) = 0.0;

            MathUtils<double>::MatrixType cof_mat = MathUtils<double>::CofactorMatrix(mat33);
            for (unsigned i = 0; i < ref33.size1(); ++i)
                for (unsigned j = 0; j < ref33.size2(); ++j)
                    KRATOS_CHECK_NEAR(cof_mat(i,j), ref33(i,j), tolerance);
        }

        /** Checks if it calculates the generalized determinant of a non-square matrix
         * Checks if it calculates the generalized determinant of a non-square matrix
         */

        KRATOS_TEST_CASE_IN_SUITE(MathUtilsGenDetMat, KratosCoreFastSuite)
        {
            constexpr double tolerance = 1e-6;

            Matrix mat23 = ZeroMatrix(2, 3);
            mat23(0,0) = 1.0;
            mat23(1,1) = 1.0;

            double det = MathUtils<double>::GeneralizedDet(mat23);

            KRATOS_CHECK_NEAR(det, 1.0, tolerance);

            Matrix mat55 = ZeroMatrix(5, 5);
            mat55(0,0) =   1.0;
            mat55(1,1) =   1.0;
            mat55(2,2) =   1.0;
            mat55(3,3) =   1.0;
            mat55(2,3) = - 1.0;
            mat55(3,2) =   1.0;
            mat55(4,4) =   2.0;

            det = MathUtils<double>::Det(mat55);

            KRATOS_CHECK_NEAR(det, 4.0, tolerance);
        }

        /** Checks if it calculates the inverse of a 1x1, 2x2, 3x3 and 4x4 matrix
         * Checks if it calculates the inverse of a 1x1, 2x2, 3x3 and 4x4 matrix
         */

        KRATOS_TEST_CASE_IN_SUITE(MathUtilsInvMat, KratosCoreFastSuite)
        {
            constexpr double tolerance = 1e-6;

            BoundedMatrix<double, 1, 1> mat11;
            mat11(0,0) = 0.896308;

            double det;
            const BoundedMatrix<double, 1, 1> inv11 = MathUtils<double>::InvertMatrix<1>(mat11, det);
            const BoundedMatrix<double, 1, 1> I11 = prod(inv11, mat11);

            KRATOS_CHECK_NEAR(I11(0,0), 1.0, tolerance);

            BoundedMatrix<double, 2, 2> mat22;
            mat22(0,0) = 0.670005;
            mat22(0,1) = 0.853367;
            mat22(1,0) = 1.47006;
            mat22(1,1) = 1.00029;

            const BoundedMatrix<double, 2, 2> inv22 = MathUtils<double>::InvertMatrix<2>(mat22, det);
            const BoundedMatrix<double, 2, 2> I22 = prod(inv22, mat22);

            for (unsigned int i = 0; i < 2; i++)
            {
                for (unsigned int j = 0; j < 2; j++)
                {
                    if (i == j)
                    {
                        KRATOS_CHECK_NEAR(I22(i,j), 1.0, tolerance);
                    }
                    else
                    {
                        KRATOS_CHECK_NEAR(I22(i,j), 0.0, tolerance);
                    }
                }
            }

            BoundedMatrix<double, 3, 3> mat33;
            mat33(0,0) = 0.678589;
            mat33(0,1) = 0.386213;
            mat33(0,2) = 0.371126;
            mat33(1,0) = 1.01524;
            mat33(1,1) = 0.403437;
            mat33(1,2) = 1.03755;
            mat33(2,0) = 0.450516;
            mat33(2,1) = 1.08225;
            mat33(2,2) = 0.972831;

            const BoundedMatrix<double, 3, 3> inv33 = MathUtils<double>::InvertMatrix<3>(mat33, det);
            const BoundedMatrix<double, 3, 3> I33 = prod(inv33, mat33);

            for (unsigned int i = 0; i < 3; i++)
            {
                for (unsigned int j = 0; j < 3; j++)
                {
                    if (i == j)
                    {
                        KRATOS_CHECK_NEAR(I33(i,j), 1.0, tolerance);
                    }
                    else
                    {
                        KRATOS_CHECK_NEAR(I33(i,j), 0.0, tolerance);
                    }
                }
            }

            BoundedMatrix<double, 4, 4> mat44;
            mat44(0,0) = 0.00959158;
            mat44(0,1) = 0.466699;
            mat44(0,2) = 0.167357;
            mat44(0,3) = 0.255465;
            mat44(1,0) = 1.6356;
            mat44(1,1) = 0.387988;
            mat44(1,2) = 1.17823;
            mat44(1,3) = 1.38661;
            mat44(2,0) = 2.57105;
            mat44(2,1) = 1.63057;
            mat44(2,2) = 2.5713;
            mat44(2,3) = 1.73297;
            mat44(3,0) = 3.40005;
            mat44(3,1) = 1.94218;
            mat44(3,2) = 2.58081;
            mat44(3,3) = 3.3083;

            const BoundedMatrix<double, 4, 4> inv44 = MathUtils<double>::InvertMatrix<4>(mat44, det);
            const BoundedMatrix<double, 4, 4> I44 = prod(inv44, mat44);

            for (unsigned int i = 0; i < 4; i++)
            {
                for (unsigned int j = 0; j < 4; j++)
                {
                    if (i == j)
                    {
                        KRATOS_CHECK_NEAR(I44(i,j), 1.0, tolerance);
                    }
                    else
                    {
                        KRATOS_CHECK_NEAR(I44(i,j), 0.0, tolerance);
                    }
                }
            }
        }

        /** Checks if it calculates the inverse of a 1x1, 2x2, 3x3 and 4x4 matrix
         * Checks if it calculates the inverse of a 1x1, 2x2, 3x3 and 4x4 matrix
         */

        KRATOS_TEST_CASE_IN_SUITE(MathUtilsInvertMatrix, KratosCoreFastSuite)
        {
            constexpr double tolerance = 1e-6;

            double det;
            Matrix inv(1,1);
            Matrix I(1,1);

            unsigned int i_dim = 1;

            Matrix mat = ZeroMatrix(i_dim, i_dim);

            mat(0,0) = 0.346432;

            MathUtils<double>::InvertMatrix(mat,inv, det);

            I = prod(inv, mat);

            for (unsigned int i = 0; i < i_dim; i++) {
                for (unsigned int j = 0; j < i_dim; j++) {
                    if (i == j) {
                        KRATOS_CHECK_NEAR(I(i,j), 1.0, tolerance);
                    } else {
                        KRATOS_CHECK_NEAR(I(i,j), 0.0, tolerance);
                    }
                }
            }

            i_dim = 2;
            mat.resize(i_dim, i_dim, false);
            inv.resize(i_dim, i_dim, false);
            I.resize(i_dim, i_dim, false);

            mat(0,0) = 0.833328;
            mat(0,1) = 0.491166;
            mat(1,0) = 0.81167;
            mat(1,1) = 1.17205;

            MathUtils<double>::InvertMatrix(mat,inv, det);

            I = prod(inv, mat);

            for (unsigned int i = 0; i < i_dim; i++) {
                for (unsigned int j = 0; j < i_dim; j++) {
                    if (i == j) {
                        KRATOS_CHECK_NEAR(I(i,j), 1.0, tolerance);
                    } else {
                        KRATOS_CHECK_NEAR(I(i,j), 0.0, tolerance);
                    }
                }
            }

            i_dim = 3;
            mat.resize(i_dim, i_dim, false);
            inv.resize(i_dim, i_dim, false);
            I.resize(i_dim, i_dim, false);

            mat(0,0) = 0.371083;
            mat(0,1) = 0.392607;
            mat(0,2) = 0.306494;
            mat(1,0) = 0.591012;
            mat(1,1) = 1.00733;
            mat(1,2) = 1.07727;
            mat(2,0) = 0.0976054;
            mat(2,1) = 2.54893;
            mat(2,2) = 1.23981;

            MathUtils<double>::InvertMatrix(mat,inv, det);

            I = prod(inv, mat);

            for (unsigned int i = 0; i < i_dim; i++) {
                for (unsigned int j = 0; j < i_dim; j++) {
                    if (i == j) {
                        KRATOS_CHECK_NEAR(I(i,j), 1.0, tolerance);
                    } else {
                        KRATOS_CHECK_NEAR(I(i,j), 0.0, tolerance);
                    }
                }
            }

            i_dim = 4;
            mat.resize(i_dim, i_dim, false);
            inv.resize(i_dim, i_dim, false);
            I.resize(i_dim, i_dim, false);

            mat(0,0) = 0.0;
            mat(0,1) = 0.979749;
            mat(0,2) = 0.494393;
            mat(0,3) = 0.23073;
            mat(1,0) = 1.79224;
            mat(1,1) = 0.198842;
            mat(1,2) = 0.074485;
            mat(1,3) = 1.45717;
            mat(2,0) = 1.6039;
            mat(2,1) = 0.673926;
            mat(2,2) = 2.63817;
            mat(2,3) = 1.0287;
            mat(3,0) = 0.366503;
            mat(3,1) = 3.02634;
            mat(3,2) = 1.24104;
            mat(3,3) = 3.62022;

            MathUtils<double>::InvertMatrix(mat,inv, det);

            I = prod(inv, mat);

            for (unsigned int i = 0; i < i_dim; i++) {
                for (unsigned int j = 0; j < i_dim; j++) {
                    if (i == j) {
                        KRATOS_CHECK_NEAR(I(i,j), 1.0, tolerance);
                    } else {
                        KRATOS_CHECK_NEAR(I(i,j), 0.0, tolerance);
                    }
                }
            }

            i_dim = 5;
            mat.resize(i_dim, i_dim, false);
            inv.resize(i_dim, i_dim, false);
            I.resize(i_dim, i_dim, false);

            mat = ZeroMatrix(5, 5);
            mat(0,0) =   1.0;
            mat(1,1) =   1.0;
            mat(2,2) =   1.0;
            mat(3,3) =   1.0;
            mat(2,3) = - 1.0;
            mat(3,2) =   1.0;
            mat(4,4) =   2.0;

            MathUtils<double>::InvertMatrix(mat,inv, det);

            KRATOS_CHECK_NEAR(det, 4.0, tolerance);

            I = prod(inv, mat);

            for (unsigned int i = 0; i < i_dim; i++) {
                for (unsigned int j = 0; j < i_dim; j++) {
                    if (i == j) {
                        KRATOS_CHECK_NEAR(I(i,j), 1.0, tolerance);
                    } else {
                        KRATOS_CHECK_NEAR(I(i,j), 0.0, tolerance);
                    }
                }
            }
        }

        /** Checks if it can solve a dense system of equations
         * Checks if it can solve a dense system of equations
         */

        KRATOS_TEST_CASE_IN_SUITE(MathUtilsSolve, KratosCoreFastSuite)
        {
            constexpr double tolerance = 1e-6;

            const std::size_t i_dim = 4;
            double det;
            Matrix A(i_dim, i_dim);
            Matrix inv(i_dim, i_dim);
            Vector b(i_dim);

            A(0,0) = 0.0;
            A(0,1) = 0.979749;
            A(0,2) = 0.494393;
            A(0,3) = 0.23073;
            A(1,0) = 1.79224;
            A(1,1) = 0.198842;
            A(1,2) = 0.074485;
            A(1,3) = 1.45717;
            A(2,0) = 1.6039;
            A(2,1) = 0.673926;
            A(2,2) = 2.63817;
            A(2,3) = 1.0287;
            A(3,0) = 0.366503;
            A(3,1) = 3.02634;
            A(3,2) = 1.24104;
            A(3,3) = 3.62022;

            b[0] = 0.0;
            b[1] = 1.0;
            b[2] = 2.0;
            b[3] = 3.0;

            MathUtils<double>::InvertMatrix(A,inv, det);

            const Vector ref_x = prod(inv, b);
            Vector x;

            MathUtils<double>::Solve(A,x,b);

            for (std::size_t i = 0; i < i_dim; i++) {
                KRATOS_CHECK_NEAR(ref_x[i], x[i], tolerance);
            }
        }

        /** Checks if it calculates correctly the inverse of a non square matrix
         * Checks if it calculates correctly the inverse of a non square matrix
         */

        KRATOS_TEST_CASE_IN_SUITE(MathUtilsGeneralizedInvertMatrix, KratosCoreFastSuite)
        {
            constexpr double tolerance = 1e-6;

            // We check the Left inverse

            const unsigned int i_dim = 2;
            const unsigned int j_dim = 3;

            Matrix mat = ZeroMatrix(i_dim, j_dim);

            mat(0,0) = 0.770724;
            mat(1,0) = 0.573294;
            mat(0,1) = 1.27699;
            mat(1,1) = 1.57776;
            mat(0,2) = 1.30216;
            mat(1,2) = 2.66483;

            double det;
            Matrix inv;

            MathUtils<double>::GeneralizedInvertMatrix(mat,inv, det);

            Matrix I = prod(mat, inv);

            for (unsigned int i = 0; i < i_dim; i++)
            {
                for (unsigned int j = 0; j < i_dim; j++)
                {
                    if (i == j)
                    {
                        KRATOS_CHECK_NEAR(I(i,j), 1.0, tolerance);
                    }
                    else
                    {
                        KRATOS_CHECK_NEAR(I(i,j), 0.0, tolerance);
                    }
                }
            }

            // We check the Right inverse
            mat.resize(j_dim, i_dim, false);
            mat = ZeroMatrix(j_dim, i_dim);

            mat(0,0) = 0.786075;
            mat(1,0) = 0.91272;
            mat(2,0) = 0.745604;
            mat(0,1) = 0.992728;
            mat(1,1) = 1.82324;
            mat(2,1) = 0.19581;

            MathUtils<double>::GeneralizedInvertMatrix(mat,inv, det);

            I = prod(inv, mat);

            for (unsigned int i = 0; i < i_dim; i++)
            {
                for (unsigned int j = 0; j < i_dim; j++)
                {
                    if (i == j)
                    {
                        KRATOS_CHECK_NEAR(I(i,j), 1.0, tolerance);
                    }
                    else
                    {
                        KRATOS_CHECK_NEAR(I(i,j), 0.0, tolerance);
                    }
                }
            }
        }

        /** Checks if it calculates the sign function
         * Checks if it calculates the sign function
         */

        KRATOS_TEST_CASE_IN_SUITE(MathUtilsSign, KratosCoreFastSuite)
        {
            int sign = MathUtils<double>::Sign(-1.0);

            KRATOS_CHECK_EQUAL(sign, -1);

            sign = MathUtils<double>::Sign(1.0);

            KRATOS_CHECK_EQUAL(sign, 1);
        }

        /** Checks if it calculates the eigen decomposition of a 3x3 system
         * Checks if it calculates the eigen decomposition of a 3x3 system
         */

        KRATOS_TEST_CASE_IN_SUITE(MathUtilsEigen, KratosCoreFastSuite)
        {
            constexpr double tolerance = 1e-6;

            BoundedMatrix<double, 3, 3> mat33;
            BoundedMatrix<double, 3, 3> eigenmat33;
            BoundedMatrix<double, 3, 3> vectormat33;

            mat33(0,0) = 0.678589;
            mat33(0,1) = 0.386213;
            mat33(0,2) = 0.371126;
            mat33(1,0) = mat33(0,1);
            mat33(1,1) = 0.403437;
            mat33(1,2) = 1.03755;
            mat33(2,0) = mat33(0,2);
            mat33(2,1) = mat33(1,2);
            mat33(2,2) = 0.972831;

            bool converged = MathUtils<double>::EigenSystem<3>(mat33, vectormat33, eigenmat33);

            BoundedMatrix<double, 3, 3> othermat33 = prod(trans(vectormat33), eigenmat33);
            BoundedMatrix<double, 3, 3> auxmat33 = prod(othermat33, vectormat33);

            for (unsigned int i = 0; i < 3; i++)
            {
                for (unsigned int j = i; j < 3; j++)
                {
                    KRATOS_CHECK_NEAR(auxmat33(i,j), mat33(i,j), tolerance);
                }
            }

            KRATOS_CHECK_EQUAL(converged, true);
        }

        /** Checks if it calculates the dot product
         * Checks if it calculates the dot product
         */

        KRATOS_TEST_CASE_IN_SUITE(MathUtilsDot, KratosCoreFastSuite)
        {
            Vector a = ZeroVector(3);
            a[1] = 1.0;
            Vector b = ZeroVector(3);
            b[0] = 1.0;

            const double c = MathUtils<double>::Dot3(a, b);
            const double d = MathUtils<double>::Dot(a, b);

            KRATOS_CHECK_EQUAL(c, 0.0);
            KRATOS_CHECK_EQUAL(d, 0.0);
        }

        /** Checks if it calculates the norm of a vector
         * Checks if it calculates the norm of a vector
         */

        KRATOS_TEST_CASE_IN_SUITE(MathUtilsNorm, KratosCoreFastSuite)
        {
            array_1d<double, 3> a = ZeroVector(3);
            a[0] = 1.0;

            const double b = MathUtils<double>::Norm3(a);
            const double c = MathUtils<double>::Norm(a);

            KRATOS_CHECK_EQUAL(b, 1.0);
            KRATOS_CHECK_EQUAL(c, 1.0);
        }

        /** Checks if it calculates the norm of a vector without underflow
         * Checks if it calculates the norm of a vector without underflow
         */

        KRATOS_TEST_CASE_IN_SUITE(MathUtilsStableNormUnderflow, KratosCoreFastSuite)
        {
            array_1d<double, 3> a = ZeroVector(3);
            a[0] = 1e-162;

            const double b = MathUtils<double>::StableNorm(a);

            KRATOS_CHECK_EQUAL(b, 1e-162);
        }

        /** Checks if it calculates the norm of a vector without overflow
         * Checks if it calculates the norm of a vector without overflow
         */

        KRATOS_TEST_CASE_IN_SUITE(MathUtilsStableNormOverflow, KratosCoreFastSuite)
        {
            array_1d<double, 3> a = ZeroVector(3);
            a[0] = 1e155;

            const double b = MathUtils<double>::StableNorm(a);

            KRATOS_CHECK_EQUAL(b, 1e155);
        }

        /** Checks if it calculates the cross product (I)
         * Checks if it calculates the cross product (I)
         */

        KRATOS_TEST_CASE_IN_SUITE(MathUtilsVectorAngleTest1, KratosCoreFastSuite)
        {
            array_1d<double, 3> a = ZeroVector(3);
            array_1d<double, 3> b = ZeroVector(3);
            a[0] = 1.0;
            b[1] = 1.0;

            const double angle = MathUtils<double>::VectorsAngle(b, a);

            KRATOS_CHECK_EQUAL(angle, Globals::Pi/2.0);
        }

        /** Checks if it calculates the cross product (II)
         * Checks if it calculates the cross product (II)
         */

        KRATOS_TEST_CASE_IN_SUITE(MathUtilsVectorAngleTest2, KratosCoreFastSuite)
        {
            array_1d<double, 3> a = ZeroVector(3);
            array_1d<double, 3> b = ZeroVector(3);
            a[0] = 1.0;
            b[0] = -1.0;

            const double angle = MathUtils<double>::VectorsAngle(b, a);

            KRATOS_CHECK_EQUAL(angle, Globals::Pi);
        }

        /** Checks if it calculates the cross product (III)
         * Checks if it calculates the cross product (III)
         */

        KRATOS_TEST_CASE_IN_SUITE(MathUtilsVectorAngleTest3, KratosCoreFastSuite)
        {
            array_1d<double, 3> a = ZeroVector(3);
            array_1d<double, 3> b = ZeroVector(3);
            a[0] = 1.0;
            a[2] = 1.0;
            a /= norm_2(a);
            b[1] = -1.0;
            b[2] = 1.0;
            b /= norm_2(b);

            const double angle = MathUtils<double>::VectorsAngle(b, a);

            KRATOS_CHECK_EQUAL(angle, Globals::Pi/3.0);
        }

        /** Checks if it calculates the angle between two vectors
         * Checks if it calculates the angle between two vectors
         */

        KRATOS_TEST_CASE_IN_SUITE(MathUtilsCross, KratosCoreFastSuite)
        {
            array_1d<double, 3> a = ZeroVector(3);
            a[1] = 2.0;
            array_1d<double, 3> b = ZeroVector(3);
            b[0] = 1.0;

            array_1d<double, 3>  c, d;

            MathUtils<double>::CrossProduct(c, b, a);
            MathUtils<double>::UnitCrossProduct(d, b, a);

            KRATOS_CHECK_EQUAL(c[2], 2.0);
            KRATOS_CHECK_EQUAL(d[2], 1.0);
        }

        /** Checks if it calculates the orthonormal base
         * Checks if it calculates the orthonormal base
         */

        KRATOS_TEST_CASE_IN_SUITE(MathUtilsOrthonormalBasis, KratosCoreFastSuite)
        {
            array_1d<double, 3> a = ZeroVector(3);
            a[1] = 1.0;

            array_1d<double, 3>  b, c;

            MathUtils<double>::OrthonormalBasisHughesMoeller(a, b, c);

            KRATOS_CHECK_EQUAL(b[0], 1.0);
            KRATOS_CHECK_EQUAL(c[2], -1.0);

            MathUtils<double>::OrthonormalBasisFrisvad(a, b, c);

            KRATOS_CHECK_EQUAL(b[0], 1.0);
            KRATOS_CHECK_EQUAL(c[2], -1.0);

            MathUtils<double>::OrthonormalBasisNaive(a, b, c);


            KRATOS_CHECK_EQUAL(b[0], 1.0);
            KRATOS_CHECK_EQUAL(c[2], -1.0);
        }

        /** Checks if it calculates the tensor product
         * Checks if it calculates the tensor product
         */

        KRATOS_TEST_CASE_IN_SUITE(MathUtilsTensor, KratosCoreFastSuite)
        {
            Vector a = ZeroVector(3);
            a[1] = 2.0;
            Vector b = ZeroVector(3);
            b[0] = 1.0;

            const Matrix c = MathUtils<double>::TensorProduct3(a, b);

            KRATOS_CHECK_EQUAL(c(0,0), 0.0);
            KRATOS_CHECK_EQUAL(c(1,0), 2.0);
            KRATOS_CHECK_EQUAL(c(0,1), 0.0);
            KRATOS_CHECK_EQUAL(c(1,1), 0.0);
        }

        /** Checks if it calculates the  matrix operations
         * Checks if it calculates the  matrix operations
         */

        KRATOS_TEST_CASE_IN_SUITE(MathUtilsMatrixOperations, KratosCoreFastSuite)
        {
            Matrix a = IdentityMatrix(3);
            Matrix b = IdentityMatrix(3);

            MathUtils<double>::AddMatrix(a, b, 0 ,0);

            KRATOS_CHECK_EQUAL(a(0,0), 2.0);
            KRATOS_CHECK_EQUAL(a(1,0), 0.0);
            KRATOS_CHECK_EQUAL(a(0,1), 0.0);
            KRATOS_CHECK_EQUAL(a(1,1), 2.0);

            MathUtils<double>::SubtractMatrix(a, b, 0 ,0);

            KRATOS_CHECK_EQUAL(a(0,0), 1.0);
            KRATOS_CHECK_EQUAL(a(1,0), 0.0);
            KRATOS_CHECK_EQUAL(a(0,1), 0.0);
            KRATOS_CHECK_EQUAL(a(1,1), 1.0);

            MathUtils<double>::WriteMatrix(a, b, 0 ,0);

            KRATOS_CHECK_EQUAL(a(0,0), 1.0);
            KRATOS_CHECK_EQUAL(a(1,0), 0.0);
            KRATOS_CHECK_EQUAL(a(0,1), 0.0);
            KRATOS_CHECK_EQUAL(a(1,1), 1.0);
        }

    } // namespace Testing
}  // namespace Kratos.

