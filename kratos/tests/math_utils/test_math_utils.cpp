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
        
        /** It test the distance function
         * It test the distance function
         */
        KRATOS_TEST_CASE_IN_SUITE(MathUtilsDistanceTest, KratosCoreMathUtilsFastSuite) 
        {
            // FIXME: I don't know what the Distnace function does
            KRATOS_ERROR << "Something went wrong" << std::endl;
        }
        
//         /** Checks if the area of the triangle is calculated correctly using Heron equation.
//          * Checks if the area of the triangle is calculated correctly using Heron equation.
//          */
//         KRATOS_TEST_CASE_IN_SUITE(MathUtilsHeronTest, KratosCoreMathUtilsFastSuite) 
//         {
//             const double area = MathUtils<double>::Heron<false>(std::sqrt(2.0), 1.0, 1.0);
// 
//             KRATOS_CHECK_NEAR(area, 0.5, TOLERANCE);
//         }
//         
//         /** Checks if it gives you the absolute value of a given value
//          * Checks if It gives you the absolute value of a given value
//          */
//         KRATOS_TEST_CASE_IN_SUITE(MathUtilsAbsTest, KratosCoreMathUtilsFastSuite) 
//         {
//             const double absolute = MathUtils<double>::Abs(-1.0);
// 
//             KRATOS_CHECK_EQUAL(absolute, 1.0);
//         }
//         
//         /** Checks if it gives you the minimum value of a given value
//          * Checks if It gives you the minimum value of a given value
//          */
//         KRATOS_TEST_CASE_IN_SUITE(MathUtilsMinTest, KratosCoreMathUtilsFastSuite) 
//         {
//             const double min = MathUtils<double>::Min(0.0,1.0);
// 
//             KRATOS_CHECK_EQUAL(min, 0.0);
//         }
//         
//         /** Checks if it gives you the maximum value of a given value
//          * Checks if It gives you the maximum value of a given value
//          */
//         KRATOS_TEST_CASE_IN_SUITE(MathUtilsMaxTest, KratosCoreMathUtilsFastSuite) 
//         {
//             const double max = MathUtils<double>::Max(0.0,1.0);
// 
//             KRATOS_CHECK_EQUAL(max, 1.0);
//         }
//         
//         /** Checks if it calculates the determinant of a 1x1, 2x2, 3x3 and 4x4 matrix 
//          * Checks if it calculates the determinant of a 1x1, 2x2, 3x3 and 4x4 matrix 
//          */
//         KRATOS_TEST_CASE_IN_SUITE(MathUtilsDetMatTest, KratosCoreMathUtilsFastSuite) 
//         {
//             boost::numeric::ublas::bounded_matrix<double, 2, 2> mat22 = ZeroMatrix(2, 2);
//             mat22(0,0) = 1.0;
//             mat22(1,1) = 1.0;
//             
//             double det = MathUtils<double>::DetMat<2>(mat22);
// 
//             KRATOS_CHECK_NEAR(det, 1.0, TOLERANCE);
//             
//             boost::numeric::ublas::bounded_matrix<double, 3, 3> mat33 = ZeroMatrix(3, 3);
//             mat33(0,0) = 1.0;
//             mat33(1,1) = 1.0;
//             mat33(2,2) = 1.0;
//             
//             det = MathUtils<double>::DetMat<3>(mat33);
// 
//             KRATOS_CHECK_NEAR(det, 1.0, TOLERANCE);
//             
//             boost::numeric::ublas::bounded_matrix<double, 4, 4> mat44 = ZeroMatrix(4, 4);
//             mat44(0,0) = 1.0;
//             mat44(1,1) = 1.0;
//             mat44(2,2) = 1.0;
//             mat44(3,3) = 1.0;
//             
//             det = MathUtils<double>::DetMat<4>(mat44);
// 
//             KRATOS_CHECK_NEAR(det, 1.0, TOLERANCE);
//         }
        
        
    } // namespace Testing
}  // namespace Kratos.

