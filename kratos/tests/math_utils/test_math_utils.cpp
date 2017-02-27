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

// External includes


// Project includes
#include "testing/testing.h"

// Utility includes
#include "utilities/math_utils.h"

namespace Kratos 
{
    namespace Testing 
    {
        constexpr double EPSILON = std::numeric_limits<double>::epsilon();
        constexpr double TOLERANCE = 1e-6;
        
        /// Tests
        
//         /** It test the distance function
//          * It test the distance function
//          */
//         
//         KRATOS_TEST_CASE_IN_SUITE(MathUtilsDistanceTest, KratosCoreMathUtilsFastSuite) 
//         {
//             // FIXME: I don't know what the Distnace function does
//         }
//         
        /** Checks if the area of the triangle is calculated correctly using Heron equation.
         * Checks if the area of the triangle is calculated correctly using Heron equation.
         */
        
        KRATOS_TEST_CASE_IN_SUITE(MathUtilsHeronTest, KratosCoreMathUtilsFastSuite) 
        {
            const double area = MathUtils<double>::Heron<false>(std::sqrt(2.0), 1.0, 1.0);

            KRATOS_CHECK_NEAR(area, 0.5, TOLERANCE);
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
            boost::numeric::ublas::bounded_matrix<double, 1, 1> mat11 = ZeroMatrix(1, 1);
            mat11(0,0) = 1.0;
            
            double det = MathUtils<double>::DetMat<1>(mat11);

            KRATOS_CHECK_NEAR(det, 1.0, TOLERANCE);
            
            boost::numeric::ublas::bounded_matrix<double, 2, 2> mat22 = ZeroMatrix(2, 2);
            mat22(0,0) = 1.0;
            mat22(1,1) = 1.0;
            
            det = MathUtils<double>::DetMat<2>(mat22);

            KRATOS_CHECK_NEAR(det, 1.0, TOLERANCE);
            
            boost::numeric::ublas::bounded_matrix<double, 3, 3> mat33 = ZeroMatrix(3, 3);
            mat33(0,0) = 1.0;
            mat33(1,1) = 1.0;
            mat33(2,2) = 1.0;
            
            det = MathUtils<double>::DetMat<3>(mat33);

            KRATOS_CHECK_NEAR(det, 1.0, TOLERANCE);
            
            boost::numeric::ublas::bounded_matrix<double, 4, 4> mat44 = ZeroMatrix(4, 4);
            mat44(0,0) = 1.0;
            mat44(1,1) = 1.0;
            mat44(2,2) = 1.0;
            mat44(3,3) = 1.0;
            
            det = MathUtils<double>::DetMat<4>(mat44);

            KRATOS_CHECK_NEAR(det, 1.0, TOLERANCE);
        }
        
        /** Checks if it calculates the inverse of a 1x1, 2x2, 3x3 and 4x4 matrix 
         * Checks if it calculates the inverse of a 1x1, 2x2, 3x3 and 4x4 matrix 
         */
        
        KRATOS_TEST_CASE_IN_SUITE(MathUtilsInvMatTest, KratosCoreMathUtilsFastSuite) 
        {
            boost::numeric::ublas::bounded_matrix<double, 1, 1> mat11 = ZeroMatrix(1, 1);
            mat11(0,0) = 2.0;
            
            double det;
            boost::numeric::ublas::bounded_matrix<double, 1, 1> inv11 = MathUtils<double>::InvertMatrix<1>(mat11, det);

            KRATOS_CHECK_NEAR(inv11(0,0), 0.5, TOLERANCE);
            
            boost::numeric::ublas::bounded_matrix<double, 2, 2> mat22 = ZeroMatrix(2, 2);
            mat22(0,0) = 1.0;
            mat22(0,1) = 2.0;
            mat22(1,0) = 3.0;
            mat22(1,1) = 4.0;
            
            boost::numeric::ublas::bounded_matrix<double, 2, 2> inv22 = MathUtils<double>::InvertMatrix<2>(mat22, det);

            KRATOS_CHECK_NEAR(inv22(0,0), -2.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv22(0,1),  1.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv22(1,0),  1.5, TOLERANCE);
            KRATOS_CHECK_NEAR(inv22(1,1), -0.5, TOLERANCE);
            
            boost::numeric::ublas::bounded_matrix<double, 3, 3> mat33 = ZeroMatrix(3, 3);
            mat33(0,0) = 1.0;
            mat33(0,1) = 2.0;
            mat33(0,2) = 2.0;
            mat33(1,0) = 2.0;
            mat33(1,1) = 2.0;
            mat33(1,2) = 2.0;
            mat33(2,0) = 2.0;
            mat33(2,1) = 2.0;
            mat33(2,2) = 3.0;
            
            boost::numeric::ublas::bounded_matrix<double, 3, 3> inv33 = MathUtils<double>::InvertMatrix<3>(mat33, det);

            KRATOS_CHECK_NEAR(inv33(0,0), -1.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv33(0,1),  1.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv33(0,2),  0.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv33(1,0),  1.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv33(1,1),  0.5, TOLERANCE);
            KRATOS_CHECK_NEAR(inv33(1,2), -1.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv33(2,0),  0.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv33(2,1), -1.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv33(2,2),  1.0, TOLERANCE);
            
            boost::numeric::ublas::bounded_matrix<double, 4, 4> mat44 = ZeroMatrix(4, 4);
            mat44(0,0) = 1.0;
            mat44(0,1) = 0.0;
            mat44(0,2) = 0.0;
            mat44(0,3) = 0.0;
            mat44(1,0) = 0.0;
            mat44(1,1) = 1.0;
            mat44(1,2) = 1.0;
            mat44(1,3) = 1.0;
            mat44(2,0) = 0.0;
            mat44(2,1) = 0.0;
            mat44(2,2) = 2.0;
            mat44(2,3) = 2.0;
            mat44(3,0) = 0.0;
            mat44(3,1) = 0.0;
            mat44(3,2) = 0.0;
            mat44(3,3) = 3.0;
            
            boost::numeric::ublas::bounded_matrix<double, 4, 4> inv44 = MathUtils<double>::InvertMatrix<4>(mat44, det);

            KRATOS_CHECK_NEAR(inv44(0,0),  1.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv44(0,1),  0.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv44(0,2),  0.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv44(0,3),  0.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv44(1,0),  0.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv44(1,1),  1.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv44(1,2), -0.5, TOLERANCE);
            KRATOS_CHECK_NEAR(inv44(1,3),  0.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv44(2,0),  0.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv44(2,1),  0.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv44(2,2),  0.5, TOLERANCE);
            KRATOS_CHECK_NEAR(inv44(2,3),  (-1.0/3.0), TOLERANCE);
            KRATOS_CHECK_NEAR(inv44(3,0),  0.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv44(3,1),  0.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv44(3,2),  0.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv44(3,3),  (1.0/3.0), TOLERANCE);
        }
        
        /** Checks if it calculates the inverse of a 1x1, 2x2, 3x3 and 4x4 matrix 
         * Checks if it calculates the inverse of a 1x1, 2x2, 3x3 and 4x4 matrix 
         */
        
        KRATOS_TEST_CASE_IN_SUITE(MathUtilsInvertMatrixTest, KratosCoreMathUtilsFastSuite) 
        {
            Matrix mat11 = ZeroMatrix(1, 1);
            mat11(0,0) = 2.0;
            
            double det;
            Matrix inv11;
            
            MathUtils<double>::InvertMatrix(mat11,inv11, det);

            KRATOS_CHECK_NEAR(inv11(0,0), 0.5, TOLERANCE);
            
            Matrix mat22 = ZeroMatrix(2, 2);
            mat22(0,0) = 1.0;
            mat22(0,1) = 2.0;
            mat22(1,0) = 3.0;
            mat22(1,1) = 4.0;
            
            Matrix inv22; 
            MathUtils<double>::InvertMatrix(mat22, inv22, det);

            KRATOS_CHECK_NEAR(inv22(0,0), -2.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv22(0,1),  1.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv22(1,0),  1.5, TOLERANCE);
            KRATOS_CHECK_NEAR(inv22(1,1), -0.5, TOLERANCE);
            
            Matrix mat33 = ZeroMatrix(3, 3);
            mat33(0,0) = 1.0;
            mat33(0,1) = 2.0;
            mat33(0,2) = 2.0;
            mat33(1,0) = 2.0;
            mat33(1,1) = 2.0;
            mat33(1,2) = 2.0;
            mat33(2,0) = 2.0;
            mat33(2,1) = 2.0;
            mat33(2,2) = 3.0;
            
            Matrix inv33;
            MathUtils<double>::InvertMatrix(mat33, inv33, det);

            KRATOS_CHECK_NEAR(inv33(0,0), -1.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv33(0,1),  1.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv33(0,2),  0.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv33(1,0),  1.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv33(1,1),  0.5, TOLERANCE);
            KRATOS_CHECK_NEAR(inv33(1,2), -1.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv33(2,0),  0.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv33(2,1), -1.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv33(2,2),  1.0, TOLERANCE);
            
            Matrix mat44 = ZeroMatrix(4, 4);
            mat44(0,0) = 1.0;
            mat44(0,1) = 0.0;
            mat44(0,2) = 0.0;
            mat44(0,3) = 0.0;
            mat44(1,0) = 0.0;
            mat44(1,1) = 1.0;
            mat44(1,2) = 1.0;
            mat44(1,3) = 1.0;
            mat44(2,0) = 0.0;
            mat44(2,1) = 0.0;
            mat44(2,2) = 2.0;
            mat44(2,3) = 2.0;
            mat44(3,0) = 0.0;
            mat44(3,1) = 0.0;
            mat44(3,2) = 0.0;
            mat44(3,3) = 3.0;
            
            Matrix inv44;
            MathUtils<double>::InvertMatrix(mat44, inv44, det);

            KRATOS_CHECK_NEAR(inv44(0,0),  1.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv44(0,1),  0.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv44(0,2),  0.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv44(0,3),  0.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv44(1,0),  0.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv44(1,1),  1.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv44(1,2), -0.5, TOLERANCE);
            KRATOS_CHECK_NEAR(inv44(1,3),  0.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv44(2,0),  0.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv44(2,1),  0.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv44(2,2),  0.5, TOLERANCE);
            KRATOS_CHECK_NEAR(inv44(2,3),  (-1.0/3.0), TOLERANCE);
            KRATOS_CHECK_NEAR(inv44(3,0),  0.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv44(3,1),  0.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv44(3,2),  0.0, TOLERANCE);
            KRATOS_CHECK_NEAR(inv44(3,3),  (1.0/3.0), TOLERANCE);
        }
        
        /** Checks if it calculates the sign function 
         * Checks if it calculates the sign function 
         */
        
        KRATOS_TEST_CASE_IN_SUITE(MathUtilsSignTest, KratosCoreMathUtilsFastSuite) 
        {
            int sign = MathUtils<double>::Sign(-1.0);
            
            KRATOS_CHECK_EQUAL(sign, -1.0);
            
            sign = MathUtils<double>::Sign(1.0);
            
            KRATOS_CHECK_EQUAL(sign, 1.0);
        }
        
        /** Checks if it calculates the eigen decomposition of a 3x3 system
         * Checks if it calculates the eigen decomposition of a 3x3 system
         */
        
        KRATOS_TEST_CASE_IN_SUITE(MathUtilsEigenTest, KratosCoreMathUtilsFastSuite) 
        {
            boost::numeric::ublas::bounded_matrix<double, 3, 3> mat33 = ZeroMatrix(3, 3);
            boost::numeric::ublas::bounded_matrix<double, 3, 3> eigenmat33;
            boost::numeric::ublas::bounded_matrix<double, 3, 3> vectormat33;
            mat33(0,0) = 1.0;
            mat33(0,1) = 0.0;
            mat33(0,2) = 1.0;
            mat33(1,0) = 0.0;
            mat33(1,1) = 2.0;
            mat33(1,2) = 0.0;
            mat33(2,0) = 1.0;
            mat33(2,1) = 0.0;
            mat33(2,2) = 1.0;
            
            bool converged = MathUtils<double>::EigenSystem<3>(mat33, vectormat33, eigenmat33);

            KRATOS_CHECK_NEAR(eigenmat33(0,0), 0.0, TOLERANCE);
            KRATOS_CHECK_NEAR(eigenmat33(1,1), 2.0, TOLERANCE);
            KRATOS_CHECK_NEAR(eigenmat33(2,2), 2.0, TOLERANCE);
            
            const double sqrt2o2 = std::sqrt(2.0)/2.0;
            
            KRATOS_CHECK_NEAR(vectormat33(0,0), sqrt2o2, TOLERANCE);
            KRATOS_CHECK_NEAR(vectormat33(0,1),  0.0, TOLERANCE);
            KRATOS_CHECK_NEAR(vectormat33(0,2), -sqrt2o2, TOLERANCE);  

            KRATOS_CHECK_NEAR(vectormat33(1,0),  0.0, TOLERANCE);
            KRATOS_CHECK_NEAR(vectormat33(1,1),  1.0, TOLERANCE);
            KRATOS_CHECK_NEAR(vectormat33(1,2),  0.0, TOLERANCE);
 
            KRATOS_CHECK_NEAR(vectormat33(2,0),  sqrt2o2, TOLERANCE);
            KRATOS_CHECK_NEAR(vectormat33(2,1),  0.0, TOLERANCE);
            KRATOS_CHECK_NEAR(vectormat33(2,2),  sqrt2o2, TOLERANCE);
            
            KRATOS_CHECK_EQUAL(converged, true);
            
        }
        
    } // namespace Testing
}  // namespace Kratos.

