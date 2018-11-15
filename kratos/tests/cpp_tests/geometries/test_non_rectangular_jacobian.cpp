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
#include "tests/cpp_tests/geometries/test_geometry.h"
#include "geometries/geometry.h"
#include "geometries/line_2d_2.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/quadrilateral_3d_4.h"

// Utility includes
#include "utilities/math_utils.h"

namespace Kratos 
{
    namespace Testing 
    {
        ///@}
        ///@name Type Definitions
        ///@{
        
        typedef GeometryData::IntegrationMethod IntegrationMethod;
        
        /// Tests
       
        /** Checks if the area of the triangle is calculated correctly using Heron equation.
         * Checks if the area of the triangle is calculated correctly using Heron equation.
         */
        
        KRATOS_TEST_CASE_IN_SUITE(LineJacobianTest, KratosNonRectangularJacobianFastSuite) 
        {
            Node<3>::Pointer PointA = GeneratePoint<Node<3>>();
            Node<3>::Pointer PointB = GeneratePoint<Node<3>>();
            
            auto geom = Line2D2<Node<3>>(PointA,PointB);

            IntegrationMethod ThisMethod = geom.GetDefaultIntegrationMethod();
            
            Matrix jacobian ( 2, 1 );
            
            geom.Jacobian( jacobian, 0, ThisMethod);
            
            const double detJ = geom.DeterminantOfJacobian(0, ThisMethod);
            
            KRATOS_CHECK_NEAR(detJ, MathUtils<double>::GeneralizedDet(jacobian), TOLERANCE);
        }
        
        /** Checks if it gives you the absolute value of a given value
         * Checks if It gives you the absolute value of a given value
         */
        
        KRATOS_TEST_CASE_IN_SUITE(TriangleJacobianTest, KratosNonRectangularJacobianFastSuite) 
        {
            Node<3>::Pointer PointA = GeneratePoint<Node<3>>();
            Node<3>::Pointer PointB = GeneratePoint<Node<3>>();
            Node<3>::Pointer PointC = GeneratePoint<Node<3>>();
            
            auto geom = Triangle3D3<Node<3>>(PointA,PointB,PointC);

            IntegrationMethod ThisMethod = geom.GetDefaultIntegrationMethod();
            
            Matrix jacobian ( 3, 2 );
            
            geom.Jacobian( jacobian, 0, ThisMethod);
            
            const double detJ = geom.DeterminantOfJacobian(0, ThisMethod);
            
            KRATOS_CHECK_NEAR(detJ, MathUtils<double>::GeneralizedDet(jacobian), TOLERANCE);
        }
        
        /** Checks if it gives you the minimum value of a given value
         * Checks if It gives you the minimum value of a given value
         */
        
        KRATOS_TEST_CASE_IN_SUITE(QuadrilateralJacobianTest, KratosNonRectangularJacobianFastSuite) 
        {
            Node<3>::Pointer PointA = Node<3>::Pointer(new Node<3>(1, 0.1, 0.2, 0.3));
            Node<3>::Pointer PointB = Node<3>::Pointer(new Node<3>(2, 0.4, 0.25, 0.35));
            Node<3>::Pointer PointC = Node<3>::Pointer(new Node<3>(3, 0.4, 0.55, 0.3));
            Node<3>::Pointer PointD = Node<3>::Pointer(new Node<3>(4, 0.15, 0.45, 0.45));
            
            auto geom = Quadrilateral3D4<Node<3>>(PointA,PointB,PointC,PointD);

            IntegrationMethod ThisMethod = geom.GetDefaultIntegrationMethod();
            
            Matrix jacobian ( 3, 2 );
            
            geom.Jacobian( jacobian, 0, ThisMethod);
            
            const double detJ = geom.DeterminantOfJacobian(0, ThisMethod);
            
            KRATOS_CHECK_NEAR(detJ, MathUtils<double>::GeneralizedDet(jacobian), TOLERANCE);
        }
    } // namespace Testing
}  // namespace Kratos.

