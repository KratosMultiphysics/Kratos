//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

// System includes
#include <limits>
#include <stdlib.h>
#include <time.h>

// External includes


// Project includes
#include "testing/testing.h"

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
        
        KRATOS_TEST_CASE_IN_SUITE(MathUtilsHeronTest, KratosCoreMathUtilsFastSuite) 
        {
            constexpr double tolerance = 1e-6;
            
            const double area = MathUtils<double>::Heron<false>(std::sqrt(2.0), 1.0, 1.0);

            KRATOS_CHECK_NEAR(area, 0.5, tolerance);
        }
        
        /** Checks if it gives you the absolute value of a given value
         * Checks if It gives you the absolute value of a given value
         */
        
        KRATOS_TEST_CASE_IN_SUITE(MathUtilsAbsTest, KratosCoreMathUtilsFastSuite) 
        {
            const double absolute = MathUtils<double>::Abs(-1.0);

            KRATOS_CHECK_EQUAL(absolute, 1.0);
        }
        
        /** Checks if it gives you the minimum value of a given value
         * Checks if It gives you the minimum value of a given value
         */
        
        KRATOS_TEST_CASE_IN_SUITE(MathUtilsMinTest, KratosCoreMathUtilsFastSuite) 
        {
            const double min = MathUtils<double>::Min(0.0,1.0);

            KRATOS_CHECK_EQUAL(min, 0.0);
        }
        
        /** Checks if it gives you the maximum value of a given value
         * Checks if It gives you the maximum value of a given value
         */
        
        KRATOS_TEST_CASE_IN_SUITE(MathUtilsMaxTest, KratosCoreMathUtilsFastSuite) 
        {
            const double max = MathUtils<double>::Max(0.0,1.0);

            KRATOS_CHECK_EQUAL(max, 1.0);
        }
        
        /** Checks if it calculates the determinant of a 1x1, 2x2, 3x3 and 4x4 matrix 
         * Checks if it calculates the determinant of a 1x1, 2x2, 3x3 and 4x4 matrix 
         */
        
        KRATOS_TEST_CASE_IN_SUITE(MathUtilsDetMatTest, KratosCoreMathUtilsFastSuite) 
        {
            constexpr double tolerance = 1e-6;
            
            boost::numeric::ublas::bounded_matrix<double, 1, 1> mat11 = ZeroMatrix(1, 1);
            mat11(0,0) = 1.0;
            
            double det = MathUtils<double>::DetMat<1>(mat11);

            KRATOS_CHECK_NEAR(det, 1.0, tolerance);
            
            boost::numeric::ublas::bounded_matrix<double, 2, 2> mat22 = ZeroMatrix(2, 2);
            mat22(0,0) = 1.0;
            mat22(1,1) = 1.0;
            
            det = MathUtils<double>::DetMat<2>(mat22);

            KRATOS_CHECK_NEAR(det, 1.0, tolerance);
            
            boost::numeric::ublas::bounded_matrix<double, 3, 3> mat33 = ZeroMatrix(3, 3);
            mat33(0,0) = 1.0;
            mat33(1,1) = 1.0;
            mat33(2,2) = 1.0;
            
            det = MathUtils<double>::DetMat<3>(mat33);

            KRATOS_CHECK_NEAR(det, 1.0, tolerance);
            
            boost::numeric::ublas::bounded_matrix<double, 4, 4> mat44 = ZeroMatrix(4, 4);
            mat44(0,0) = 1.0;
            mat44(1,1) = 1.0;
            mat44(2,2) = 1.0;
            mat44(3,3) = 1.0;
            
            det = MathUtils<double>::DetMat<4>(mat44);

            KRATOS_CHECK_NEAR(det, 1.0, tolerance);
        }
        
        /** Checks if it calculates the generalized determinant of a non-square matrix
         * Checks if it calculates the generalized determinant of a non-square matrix
         */
        
        KRATOS_TEST_CASE_IN_SUITE(MathUtilsGenDetMatTest, KratosCoreMathUtilsFastSuite) 
        {
            constexpr double tolerance = 1e-6;
            
            Matrix mat23 = ZeroMatrix(2, 3);
            mat23(0,0) = 1.0;
            mat23(1,1) = 1.0;
            
            double det = MathUtils<double>::GeneralizedDet(mat23);
            
            KRATOS_CHECK_NEAR(det, 0.0, tolerance);
            
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
        
        KRATOS_TEST_CASE_IN_SUITE(MathUtilsInvMatTest, KratosCoreMathUtilsFastSuite) 
        {
            constexpr double tolerance = 1e-6;
            
            srand (1);

            boost::numeric::ublas::bounded_matrix<double, 1, 1> mat11;
            mat11(0,0) = rand() % 10 + 1;
            
            double det;
            const boost::numeric::ublas::bounded_matrix<double, 1, 1> inv11 = MathUtils<double>::InvertMatrix<1>(mat11, det);
            const boost::numeric::ublas::bounded_matrix<double, 1, 1> I11 = prod(inv11, mat11);
            
            KRATOS_CHECK_NEAR(I11(0,0), 1.0, tolerance);
            
            boost::numeric::ublas::bounded_matrix<double, 2, 2> mat22;
            for (unsigned int i = 0; i < 2; i++)
            {
                for (unsigned int j = 0; j < 2; j++)
                {
                    mat22(i, j)= rand() % 10 + 1;
                }
            }
            
            const boost::numeric::ublas::bounded_matrix<double, 2, 2> inv22 = MathUtils<double>::InvertMatrix<2>(mat22, det);
            const boost::numeric::ublas::bounded_matrix<double, 2, 2> I22 = prod(inv22, mat22);
            
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
            
            boost::numeric::ublas::bounded_matrix<double, 3, 3> mat33;
            for (unsigned int i = 0; i < 3; i++)
            {
                for (unsigned int j = 0; j < 3; j++)
                {
                    mat33(i, j)= rand() % 10 + 1;
                }
            }
            
            const boost::numeric::ublas::bounded_matrix<double, 3, 3> inv33 = MathUtils<double>::InvertMatrix<3>(mat33, det);
            const boost::numeric::ublas::bounded_matrix<double, 3, 3> I33 = prod(inv33, mat33);
            
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
            
            boost::numeric::ublas::bounded_matrix<double, 4, 4> mat44;
            for (unsigned int i = 0; i < 4; i++)
            {
                for (unsigned int j = 0; j < 4; j++)
                {
                    mat44(i, j)= rand() % 10 + 1;
                }
            }
            
            const boost::numeric::ublas::bounded_matrix<double, 4, 4> inv44 = MathUtils<double>::InvertMatrix<4>(mat44, det);
            const boost::numeric::ublas::bounded_matrix<double, 4, 4> I44 = prod(inv44, mat44);
            
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
        
        KRATOS_TEST_CASE_IN_SUITE(MathUtilsInvertMatrixTest, KratosCoreMathUtilsFastSuite) 
        {
            constexpr double tolerance = 1e-6;
            
            srand (1);
            
            for (unsigned int i_dim = 1; i_dim <= 4; i_dim++)
            {
                Matrix mat = ZeroMatrix(i_dim, i_dim);
                
                for (unsigned int i = 0; i < i_dim; i++)
                {
                    for (unsigned int j = 0; j < i_dim; j++)
                    {
                        mat(i, j)= rand() % 10 + 1;
                    }
                }
                
                double det;
                Matrix inv;
                
                MathUtils<double>::InvertMatrix(mat,inv, det);
                
                const Matrix I = prod(inv, mat);
                
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
        }
        
        /** Checks if it calculates correctly the inverse of a non square matrix
         * Checks if it calculates correctly the inverse of a non square matrix
         */
        
        KRATOS_TEST_CASE_IN_SUITE(MathUtilsGeneralizedInvertMatrixTest, KratosCoreMathUtilsFastSuite) 
        {
            constexpr double tolerance = 1e-6;
            
            srand (1);
            
            // We check the Left inverse
            for (unsigned int i_dim = 1; i_dim <= 4; i_dim++)
            {
                Matrix mat = ZeroMatrix(i_dim, i_dim - 1);
                
                for (unsigned int i = 0; i < i_dim; i++)
                {
                    for (unsigned int j = 0; j < i_dim - 1; j++)
                    {
                        mat(i, j)= rand() % 10 + 1;
                    }
                }
                
                double det;
                Matrix inv;
                
                MathUtils<double>::GeneralizedInvertMatrix(mat,inv, det);
                
                const Matrix I = prod(inv, mat);
                
                for (unsigned int i = 0; i < i_dim - 1; i++)
                {
                    for (unsigned int j = 0; j < i_dim - 1; j++)
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
            
            // We check the Right
            for (unsigned int i_dim = 1; i_dim <= 4; i_dim++)
            {
                Matrix mat = ZeroMatrix(i_dim - 1, i_dim);
                
                for (unsigned int i = 0; i < i_dim - 1; i++)
                {
                    for (unsigned int j = 0; j < i_dim; j++)
                    {
                        mat(i, j)= rand() % 10 + 1;
                    }
                }
                
                double det;
                Matrix inv;
                
                MathUtils<double>::GeneralizedInvertMatrix(mat,inv, det);
                
                const Matrix I = prod(inv, mat);
                
                for (unsigned int i = 0; i < i_dim - 1; i++)
                {
                    for (unsigned int j = 0; j < i_dim - 1; j++)
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
        }
        
        /** Checks if it calculates the sign function 
         * Checks if it calculates the sign function 
         */
        
        KRATOS_TEST_CASE_IN_SUITE(MathUtilsSignTest, KratosCoreMathUtilsFastSuite) 
        {
            int sign = MathUtils<double>::Sign(-1.0);
            
            KRATOS_CHECK_EQUAL(sign, -1);
            
            sign = MathUtils<double>::Sign(1.0);
            
            KRATOS_CHECK_EQUAL(sign, 1);
        }
        
        /** Checks if it calculates the eigen decomposition of a 3x3 system
         * Checks if it calculates the eigen decomposition of a 3x3 system
         */
        
        KRATOS_TEST_CASE_IN_SUITE(MathUtilsEigenTest, KratosCoreMathUtilsFastSuite) 
        {
            constexpr double tolerance = 1e-6;
            
            boost::numeric::ublas::bounded_matrix<double, 3, 3> mat33;
            boost::numeric::ublas::bounded_matrix<double, 3, 3> eigenmat33;
            boost::numeric::ublas::bounded_matrix<double, 3, 3> vectormat33;
            
            srand (1);
            
            for (unsigned int i = 0; i < 3; i++)
            {
                for (unsigned int j = i; j < 3; j++)
                {
                    mat33(i, j) = rand() % 10 + 1;
                    mat33(j, i) = mat33(i, j);
                }
            }
            
            bool converged = MathUtils<double>::EigenSystem<3>(mat33, vectormat33, eigenmat33);

            boost::numeric::ublas::bounded_matrix<double, 3, 3> auxmat33 = prod(trans(vectormat33), eigenmat33);
            auxmat33 = prod(auxmat33, vectormat33);
            
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
        
        KRATOS_TEST_CASE_IN_SUITE(MathUtilsDotTest, KratosCoreMathUtilsFastSuite) 
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
        
        KRATOS_TEST_CASE_IN_SUITE(MathUtilsNormTest, KratosCoreMathUtilsFastSuite) 
        {
            array_1d<double, 3> a = ZeroVector(3);
            a[0] = 1.0;

            const double b = MathUtils<double>::Norm3(a);
            const double c = MathUtils<double>::Norm(a);
            
            KRATOS_CHECK_EQUAL(b, 1.0);
            KRATOS_CHECK_EQUAL(c, 1.0);
        }
        
        /** Checks if it calculates the cross product 
         * Checks if it calculates the cross product 
         */
        
        KRATOS_TEST_CASE_IN_SUITE(MathUtilsCrossTest, KratosCoreMathUtilsFastSuite) 
        {
            array_1d<double, 3> a = ZeroVector(3);
            a[1] = 2.0;
            array_1d<double, 3> b = ZeroVector(3);
            b[0] = 1.0;

            const array_1d<double, 3>  c = MathUtils<double>::CrossProduct(a, b);
            const array_1d<double, 3>  d = MathUtils<double>::UnitCrossProduct(a, b);
            
            KRATOS_CHECK_EQUAL(c[2], 2.0);
            KRATOS_CHECK_EQUAL(d[2], 1.0);
        }
        
        /** Checks if it calculates the tensor product 
         * Checks if it calculates the tensor product 
         */
        
        KRATOS_TEST_CASE_IN_SUITE(MathUtilsTensorTest, KratosCoreMathUtilsFastSuite) 
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
        
        KRATOS_TEST_CASE_IN_SUITE(MathUtilsMatrixOperationsTest, KratosCoreMathUtilsFastSuite) 
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

