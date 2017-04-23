// KRATOS  ___|  |       |       |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//           | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License: BSD License
//   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Vicente Mataix Ferrándiz
//

// System includes
#include <set>

// External includes

// Project includes
#include "testing/testing.h"
#include "spaces/ublas_space.h"
#include "includes/properties.h"
#include "includes/model_part.h"
#include "utilities/math_utils.h"

/* TRIANGLES */
#include "geometries/triangle_2d_3.h"
/* QUADRILATERALS */
#include "geometries/quadrilateral_2d_4.h"

/* GAUSS-LEGENDRE */
#include "integration/triangle_gauss_legendre_integration_points.h"
#include "integration/quadrilateral_gauss_legendre_integration_points.h"

namespace Kratos 
{
    namespace Testing 
    {

        typedef Point<3>                                             PointType;
        typedef Node<3>                                               NodeType;
        typedef Geometry<NodeType>                            GeometryNodeType;
        typedef Geometry<PointType>                          GeometryPointType;
        
        ///Type definition for integration methods
        typedef GeometryData::IntegrationMethod              IntegrationMethod;
        typedef IntegrationPoint<2>                       IntegrationPointType;
        typedef GeometryNodeType::IntegrationPointsArrayType IntegrationPointsType;

        /** 
         * Checks if the criteria for computing the integral is the correct one
         * Checks mass matrix computed
         */
    
        KRATOS_TEST_CASE_IN_SUITE(TestMassMatrixIntegrationTriangle, ContactStructuralApplicationFastSuite)
        {
            ModelPart ModelPart("Main");
            
            // First we create the nodes 
            NodeType::Pointer pNode1 = ModelPart.CreateNewNode(0,-0.2,0.1,0.0);
            NodeType::Pointer pNode2 = ModelPart.CreateNewNode(1,1.0,0.1,0.0);
            NodeType::Pointer pNode3 = ModelPart.CreateNewNode(2,0.2,1.2,0.0);
            NodeType::Pointer pNode4 = ModelPart.CreateNewNode(3,0.6,0.4,0.0);
            
            // Now we create the "conditions"
            std::vector<NodeType::Pointer> ConditionNodes0 (3);
            
            ConditionNodes0[0] = pNode1;
            ConditionNodes0[1] = pNode2;
            ConditionNodes0[2] = pNode3;
            
            Triangle2D3 <Node<3>> triangle0( ConditionNodes0 );
            
            std::vector<NodeType::Pointer> ConditionNodes1 (3);
            
            ConditionNodes1[0] = pNode1;
            ConditionNodes1[1] = pNode2;
            ConditionNodes1[2] = pNode4;
            
            Triangle2D3 <Node<3>> triangle1( ConditionNodes1 );
            
            std::vector<NodeType::Pointer> ConditionNodes2 (3);
            
            ConditionNodes2[0] = pNode2;
            ConditionNodes2[1] = pNode3;
            ConditionNodes2[2] = pNode4;
            
            Triangle2D3 <Node<3>> triangle2( ConditionNodes2 );
            
            std::vector<NodeType::Pointer> ConditionNodes3 (3);
            
            ConditionNodes3[0] = pNode3;
            ConditionNodes3[1] = pNode1;
            ConditionNodes3[2] = pNode4;
            
            Triangle2D3 <Node<3>> triangle3( ConditionNodes3 );
            
            // We calculate the integral of the mass matrix (assuming constant density)
            GeometryNodeType::IntegrationPointsArrayType IntegrationPoints = Quadrature<TriangleGaussLegendreIntegrationPoints2, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            
            boost::numeric::ublas::bounded_matrix<double, 3, 3> MassMatrix0 = ZeroMatrix(3, 3);
            
            for (unsigned int PointNumber = 0; PointNumber < IntegrationPoints.size(); PointNumber++)
            {
                Vector N;
                const PointType& LocalPoint = IntegrationPoints[PointNumber].Coordinates();
                triangle0.ShapeFunctionsValues( N, LocalPoint );
                const double DetJ = triangle0.DeterminantOfJacobian( LocalPoint );
                const double weight = IntegrationPoints[PointNumber].Weight();
                
                for (unsigned int inode = 0; inode < 3; inode++)
                {
                    for (unsigned int jnode = 0; jnode < 3; jnode++)
                    {
                        MassMatrix0(inode, jnode) += DetJ * weight * N(inode) * N(jnode);
                    }
                }
            }
     
            boost::numeric::ublas::bounded_matrix<double, 3, 3> MassMatrix1 = ZeroMatrix(3, 3);
            
            for (unsigned int PointNumber = 0; PointNumber < IntegrationPoints.size(); PointNumber++)
            {
                Vector N1;
                Vector N2;
                Vector N3;
                
                const PointType& LocalPoint0 = IntegrationPoints[PointNumber].Coordinates();
                
                PointType gp_global;
                
                triangle1.GlobalCoordinates(gp_global, LocalPoint0);
                PointType LocalPoint1;
                triangle0.PointLocalCoordinates(LocalPoint1, gp_global);
                
                triangle2.GlobalCoordinates(gp_global, LocalPoint0);
                PointType LocalPoint2; 
                triangle0.PointLocalCoordinates(LocalPoint2, gp_global);
                
                triangle3.GlobalCoordinates(gp_global, LocalPoint0);
                PointType LocalPoint3; 
                triangle0.PointLocalCoordinates(LocalPoint3, gp_global);
                
                triangle0.ShapeFunctionsValues( N1, LocalPoint1 );
                triangle0.ShapeFunctionsValues( N2, LocalPoint2 );
                triangle0.ShapeFunctionsValues( N3, LocalPoint3 );
                
                const double DetJ1  = triangle1.DeterminantOfJacobian( LocalPoint0 );
                const double DetJ2  = triangle2.DeterminantOfJacobian( LocalPoint0 );
                const double DetJ3  = triangle3.DeterminantOfJacobian( LocalPoint0 );
                
                const double weight = IntegrationPoints[PointNumber].Weight();
                
                for (unsigned int inode = 0; inode < 3; inode++)
                {
                    for (unsigned int jnode = 0; jnode < 3; jnode++)
                    {
                        MassMatrix1(inode, jnode) += DetJ1 * weight * N1[inode] * N1[jnode] \
                                                   + DetJ2 * weight * N2[inode] * N2[jnode] \
                                                   + DetJ3 * weight * N3[inode] * N3[jnode];
                    }
                }
            }
            
            const double Tolerance = 1.0e-6;
            for (unsigned int inode = 0; inode < 3; inode++)
            {
                for (unsigned int jnode = 0; jnode < 3; jnode++)
                {
                    KRATOS_CHECK_NEAR(MassMatrix0(inode,jnode), MassMatrix1(inode,jnode), Tolerance);
                }
            }
        }
        
        /** 
         * Checks if the criteria for computing the integral is the correct one
         * Checks mass matrix computed
         */
    
        KRATOS_TEST_CASE_IN_SUITE(TestMassMatrixIntegrationQuadrilateral, ContactStructuralApplicationFastSuite)
        {
            ModelPart ModelPart("Main");
            
            // First we create the nodes 
            NodeType::Pointer pNode1 = ModelPart.CreateNewNode(0,   0.0,  0.0, 0.0);
            NodeType::Pointer pNode2 = ModelPart.CreateNewNode(1,   1.0,- 0.1, 0.0);
            NodeType::Pointer pNode3 = ModelPart.CreateNewNode(2,   1.2,  1.1, 0.0);
            NodeType::Pointer pNode4 = ModelPart.CreateNewNode(3, - 0.1,  1.3, 0.0);
            
            // Now we create the "conditions"
            std::vector<NodeType::Pointer> ConditionNodes0 (4);
            
            ConditionNodes0[0] = pNode1;
            ConditionNodes0[1] = pNode2;
            ConditionNodes0[2] = pNode3;
            ConditionNodes0[3] = pNode4;
            
            Quadrilateral2D4 <Node<3>> quadrilateral0( ConditionNodes0 );
            
            std::vector<NodeType::Pointer> ConditionNodes1 (3);
            
            ConditionNodes1[0] = pNode1;
            ConditionNodes1[1] = pNode2;
            ConditionNodes1[2] = pNode3;
            
            Triangle2D3 <Node<3>> triangle1( ConditionNodes1 );
            
            std::vector<NodeType::Pointer> ConditionNodes2 (3);
            
            ConditionNodes2[0] = pNode1;
            ConditionNodes2[1] = pNode3;
            ConditionNodes2[2] = pNode4;
            
            Triangle2D3 <Node<3>> triangle2( ConditionNodes2 );
            
            // We calculate the integral of the mass matrix (assuming constant density)
            GeometryNodeType::IntegrationPointsArrayType IntegrationPointsQuadrilateral = Quadrature<QuadrilateralGaussLegendreIntegrationPoints2, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            GeometryNodeType::IntegrationPointsArrayType IntegrationPointsTriangle = Quadrature<TriangleGaussLegendreIntegrationPoints5, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            
            boost::numeric::ublas::bounded_matrix<double, 4, 4> MassMatrix0 = ZeroMatrix(4, 4);
            
            for (unsigned int PointNumber = 0; PointNumber < IntegrationPointsQuadrilateral.size(); PointNumber++)
            {
                Vector N;
                const PointType& LocalPoint = IntegrationPointsQuadrilateral[PointNumber].Coordinates();
                quadrilateral0.ShapeFunctionsValues( N, LocalPoint );
                const double DetJ = quadrilateral0.DeterminantOfJacobian( LocalPoint );
                const double weight = IntegrationPointsQuadrilateral[PointNumber].Weight();
                
                for (unsigned int inode = 0; inode < 4; inode++)
                {
                    for (unsigned int jnode = 0; jnode < 4; jnode++)
                    {
                        MassMatrix0(inode, jnode) += DetJ * weight * N[inode] * N[jnode];
                    }
                }
            }
     
            boost::numeric::ublas::bounded_matrix<double, 4, 4> MassMatrix1 = ZeroMatrix(4, 4);
            
            for (unsigned int PointNumber = 0; PointNumber < IntegrationPointsTriangle.size(); PointNumber++)
            {
                Vector N1;
                Vector N2;
                
                const PointType& LocalPoint0 = IntegrationPointsTriangle[PointNumber].Coordinates();
                
                PointType gp_global;
                
                triangle1.GlobalCoordinates(gp_global, LocalPoint0);
                PointType LocalPoint1;
                quadrilateral0.PointLocalCoordinates(LocalPoint1, gp_global);
                
                triangle2.GlobalCoordinates(gp_global, LocalPoint0);
                PointType LocalPoint2; 
                quadrilateral0.PointLocalCoordinates(LocalPoint2, gp_global);
                
                quadrilateral0.ShapeFunctionsValues( N1, LocalPoint1 );
                quadrilateral0.ShapeFunctionsValues( N2, LocalPoint2 );
                
                const double DetJ1  = triangle1.DeterminantOfJacobian( LocalPoint0 );
                const double DetJ2  = triangle2.DeterminantOfJacobian( LocalPoint0 );
                
                const double weight = IntegrationPointsTriangle[PointNumber].Weight();
                
                for (unsigned int inode = 0; inode < 4; inode++)
                {
                    for (unsigned int jnode = 0; jnode < 4; jnode++)
                    {                        
                        MassMatrix1(inode, jnode ) += DetJ1 * weight * N1[inode] * N1[jnode] 
                                                    + DetJ2 * weight * N2[inode] * N2[jnode];
                    }
                }
            }

            
            // Debug
//             KRATOS_WATCH(MassMatrix0)
//             KRATOS_WATCH(MassMatrix1)
            
            const double Tolerance = 1.0e-6;
            for (unsigned int inode = 0; inode < 4; inode++)
            {
                for (unsigned int jnode = 0; jnode < 4; jnode++)
                {
                    KRATOS_CHECK_NEAR(MassMatrix0(inode,jnode), MassMatrix1(inode,jnode), Tolerance);
                }
            }
        }
        
        /** 
         * Checks if the criteria for computing the integral is the correct one
         * Checks mass matrix computed
         */
    
        KRATOS_TEST_CASE_IN_SUITE(TestMassMatrixIntegrationQuadrilateralDeformed, ContactStructuralApplicationFastSuite)
        {
            ModelPart ModelPart("Main");
            
            // First we create the nodes 
            NodeType::Pointer pNode0 = ModelPart.CreateNewNode(0,   0.5,  0.4, 0.0);
            NodeType::Pointer pNode1 = ModelPart.CreateNewNode(1,   0.0,  0.0, 0.0);
            NodeType::Pointer pNode2 = ModelPart.CreateNewNode(2,   1.0,- 0.1, 0.0);
            NodeType::Pointer pNode3 = ModelPart.CreateNewNode(3,   1.2,  1.1, 0.0);
            NodeType::Pointer pNode4 = ModelPart.CreateNewNode(4, - 0.1,  1.3, 0.0);
            
            // Now we create the "conditions"
            std::vector<NodeType::Pointer> ConditionNodes0 (4);
            
            ConditionNodes0[0] = pNode1;
            ConditionNodes0[1] = pNode2;
            ConditionNodes0[2] = pNode3;
            ConditionNodes0[3] = pNode4;
            
            Quadrilateral2D4 <Node<3>> quadrilateral0( ConditionNodes0 );
            
            std::vector<NodeType::Pointer> ConditionNodes1 (3);
            
            ConditionNodes1[0] = pNode1;
            ConditionNodes1[1] = pNode2;
            ConditionNodes1[2] = pNode0;
            
            Triangle2D3 <Node<3>> triangle1( ConditionNodes1 );
            
            std::vector<NodeType::Pointer> ConditionNodes2 (3);
            
            ConditionNodes2[0] = pNode2;
            ConditionNodes2[1] = pNode3;
            ConditionNodes2[2] = pNode0;
            
            Triangle2D3 <Node<3>> triangle2( ConditionNodes2 );
            
            std::vector<NodeType::Pointer> ConditionNodes3 (3);
            
            ConditionNodes3[0] = pNode3;
            ConditionNodes3[1] = pNode4;
            ConditionNodes3[2] = pNode0;
            
            Triangle2D3 <Node<3>> triangle3( ConditionNodes3 );
            
            std::vector<NodeType::Pointer> ConditionNodes4 (3);
            
            ConditionNodes4[0] = pNode4;
            ConditionNodes4[1] = pNode1;
            ConditionNodes4[2] = pNode0;
            
            Triangle2D3 <Node<3>> triangle4( ConditionNodes4 );
            
            // We calculate the integral of the mass matrix (assuming constant density)
            GeometryNodeType::IntegrationPointsArrayType IntegrationPointsQuadrilateral = Quadrature<QuadrilateralGaussLegendreIntegrationPoints2, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            GeometryNodeType::IntegrationPointsArrayType IntegrationPointsTriangle = Quadrature<TriangleGaussLegendreIntegrationPoints5, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            
            boost::numeric::ublas::bounded_matrix<double, 4, 4> MassMatrix0 = ZeroMatrix(4, 4);
            
            for (unsigned int PointNumber = 0; PointNumber < IntegrationPointsQuadrilateral.size(); PointNumber++)
            {
                Vector N;
                const PointType& LocalPoint = IntegrationPointsQuadrilateral[PointNumber].Coordinates();
                quadrilateral0.ShapeFunctionsValues( N, LocalPoint );
                const double DetJ = quadrilateral0.DeterminantOfJacobian( LocalPoint );
                const double weight = IntegrationPointsQuadrilateral[PointNumber].Weight();
                
                for (unsigned int inode = 0; inode < 4; inode++)
                {
                    for (unsigned int jnode = 0; jnode < 4; jnode++)
                    {
                        MassMatrix0(inode, jnode) += DetJ * weight * N[inode] * N[jnode];
                    }
                }
            }
     
            boost::numeric::ublas::bounded_matrix<double, 4, 4> MassMatrix1 = ZeroMatrix(4, 4);
            
            for (unsigned int PointNumber = 0; PointNumber < IntegrationPointsTriangle.size(); PointNumber++)
            {
                Vector N1;
                Vector N2;
                Vector N3;
                Vector N4;
                
                const PointType& LocalPoint0 = IntegrationPointsTriangle[PointNumber].Coordinates();
                
                PointType gp_global;
                
                triangle1.GlobalCoordinates(gp_global, LocalPoint0);
                PointType LocalPoint1;
                quadrilateral0.PointLocalCoordinates(LocalPoint1, gp_global);
                
                triangle2.GlobalCoordinates(gp_global, LocalPoint0);
                PointType LocalPoint2; 
                quadrilateral0.PointLocalCoordinates(LocalPoint2, gp_global);
                
                triangle3.GlobalCoordinates(gp_global, LocalPoint0);
                PointType LocalPoint3; 
                quadrilateral0.PointLocalCoordinates(LocalPoint3, gp_global);
                
                triangle4.GlobalCoordinates(gp_global, LocalPoint0);
                PointType LocalPoint4; 
                quadrilateral0.PointLocalCoordinates(LocalPoint4, gp_global);
                
                quadrilateral0.ShapeFunctionsValues( N1, LocalPoint1 );
                quadrilateral0.ShapeFunctionsValues( N2, LocalPoint2 );
                quadrilateral0.ShapeFunctionsValues( N3, LocalPoint3 );
                quadrilateral0.ShapeFunctionsValues( N4, LocalPoint4 );
                
                const double DetJ1  = triangle1.DeterminantOfJacobian( LocalPoint0 );
                const double DetJ2  = triangle2.DeterminantOfJacobian( LocalPoint0 );
                const double DetJ3  = triangle3.DeterminantOfJacobian( LocalPoint0 );
                const double DetJ4  = triangle4.DeterminantOfJacobian( LocalPoint0 );
                
                const double weight = IntegrationPointsTriangle[PointNumber].Weight();
                
                for (unsigned int inode = 0; inode < 4; inode++)
                {
                    for (unsigned int jnode = 0; jnode < 4; jnode++)
                    {                        
                        MassMatrix1(inode, jnode ) += DetJ1 * weight * N1[inode] * N1[jnode] 
                                                   +  DetJ2 * weight * N2[inode] * N2[jnode]
                                                   +  DetJ3 * weight * N3[inode] * N3[jnode]
                                                   +  DetJ4 * weight * N4[inode] * N4[jnode];
                    }
                }
            }
                        
//             // Debug
//             KRATOS_WATCH(MassMatrix0)
//             KRATOS_WATCH(MassMatrix1)
            
            const double Tolerance = 1.0e-6;
            for (unsigned int inode = 0; inode < 4; inode++)
            {
                for (unsigned int jnode = 0; jnode < 4; jnode++)
                {
                    KRATOS_CHECK_NEAR(MassMatrix0(inode,jnode), MassMatrix1(inode,jnode), Tolerance);
                }
            }
            
            array_1d<double, 3> DispArray = ZeroVector(3);
            DispArray[0] = 0.1;
            DispArray[1] = 0.2;
            
            pNode0->Coordinates() += 0.5  * DispArray;
            pNode1->Coordinates() += 1.0  * DispArray;
            pNode2->Coordinates() += 1.5  * DispArray;
            pNode3->Coordinates() += 0.75 * DispArray;
            pNode4->Coordinates() += 1.0  * DispArray;
            
            MassMatrix0 = ZeroMatrix(4, 4);
            
            for (unsigned int PointNumber = 0; PointNumber < IntegrationPointsQuadrilateral.size(); PointNumber++)
            {
                Vector N;
                const PointType& LocalPoint = IntegrationPointsQuadrilateral[PointNumber].Coordinates();
                quadrilateral0.ShapeFunctionsValues( N, LocalPoint );
                const double DetJ = quadrilateral0.DeterminantOfJacobian( LocalPoint );
                const double weight = IntegrationPointsQuadrilateral[PointNumber].Weight();
                
                for (unsigned int inode = 0; inode < 4; inode++)
                {
                    for (unsigned int jnode = 0; jnode < 4; jnode++)
                    {
                        MassMatrix0(inode, jnode) += DetJ * weight * N[inode] * N[jnode];
                    }
                }
            }
            
            MassMatrix1 = ZeroMatrix(4, 4);
            
            for (unsigned int PointNumber = 0; PointNumber < IntegrationPointsTriangle.size(); PointNumber++)
            {
                Vector N1;
                Vector N2;
                Vector N3;
                Vector N4;
                
                const PointType& LocalPoint0 = IntegrationPointsTriangle[PointNumber].Coordinates();
                
                PointType gp_global;
                
                triangle1.GlobalCoordinates(gp_global, LocalPoint0);
                PointType LocalPoint1;
                quadrilateral0.PointLocalCoordinates(LocalPoint1, gp_global);
                
                triangle2.GlobalCoordinates(gp_global, LocalPoint0);
                PointType LocalPoint2; 
                quadrilateral0.PointLocalCoordinates(LocalPoint2, gp_global);
                
                triangle3.GlobalCoordinates(gp_global, LocalPoint0);
                PointType LocalPoint3; 
                quadrilateral0.PointLocalCoordinates(LocalPoint3, gp_global);
                
                triangle4.GlobalCoordinates(gp_global, LocalPoint0);
                PointType LocalPoint4; 
                quadrilateral0.PointLocalCoordinates(LocalPoint4, gp_global);
                
                quadrilateral0.ShapeFunctionsValues( N1, LocalPoint1 );
                quadrilateral0.ShapeFunctionsValues( N2, LocalPoint2 );
                quadrilateral0.ShapeFunctionsValues( N3, LocalPoint3 );
                quadrilateral0.ShapeFunctionsValues( N4, LocalPoint4 );
                
                const double DetJ1 = triangle1.DeterminantOfJacobian( LocalPoint0 );
                const double DetJ2 = triangle2.DeterminantOfJacobian( LocalPoint0 );
                const double DetJ3 = triangle3.DeterminantOfJacobian( LocalPoint0 );
                const double DetJ4 = triangle4.DeterminantOfJacobian( LocalPoint0 );
                
                const double weight = IntegrationPointsTriangle[PointNumber].Weight();
                
                for (unsigned int inode = 0; inode < 4; inode++)
                {
                    for (unsigned int jnode = 0; jnode < 4; jnode++)
                    {                        
                        MassMatrix1(inode, jnode ) += DetJ1 * weight * N1[inode] * N1[jnode] 
                                                   +  DetJ2 * weight * N2[inode] * N2[jnode]
                                                   +  DetJ3 * weight * N3[inode] * N3[jnode]
                                                   +  DetJ4 * weight * N4[inode] * N4[jnode];
                    }
                }
            }
                        
//             // Debug
//             KRATOS_WATCH(MassMatrix0)
//             KRATOS_WATCH(MassMatrix1)
            
            for (unsigned int inode = 0; inode < 4; inode++)
            {
                for (unsigned int jnode = 0; jnode < 4; jnode++)
                {
                    KRATOS_CHECK_NEAR(MassMatrix0(inode,jnode), MassMatrix1(inode,jnode), Tolerance);
                }
            }
            
        }
        
//         /** NOTE: Doesn't work, which means that this method produces error.
//          * Checks if the criteria for computing the integral is the correct one interpolating the local coordinates
//          * Checks mass matrix computed
//          */
//         
//         KRATOS_TEST_CASE_IN_SUITE(TestMassMatrixIntegrationQuadrilateralInterpolated, ContactStructuralApplicationFastSuite)
//         {
//             ModelPart ModelPart("Main");
//             
//             // First we create the nodes 
//             NodeType::Pointer pNode1 = ModelPart.CreateNewNode(0,   0.0,  0.0, 0.0);
//             NodeType::Pointer pNode2 = ModelPart.CreateNewNode(1,   1.0,- 0.1, 0.0);
//             NodeType::Pointer pNode3 = ModelPart.CreateNewNode(2,   1.2,  1.1, 0.0);
//             NodeType::Pointer pNode4 = ModelPart.CreateNewNode(3, - 0.1,  1.3, 0.0);
//             
//             // Now we create the "conditions"
//             std::vector<NodeType::Pointer> ConditionNodes0 (4);
//             
//             ConditionNodes0[0] = pNode1;
//             ConditionNodes0[1] = pNode2;
//             ConditionNodes0[2] = pNode3;
//             ConditionNodes0[3] = pNode4;
//             
//             Quadrilateral2D4 <Node<3>> quadrilateral0( ConditionNodes0 );
//             
//             std::vector<NodeType::Pointer> ConditionNodes1 (3);
//             
//             ConditionNodes1[0] = pNode1;
//             ConditionNodes1[1] = pNode2;
//             ConditionNodes1[2] = pNode3;
//             
//             Triangle2D3 <Node<3>> triangle1( ConditionNodes1 );
//             
//             std::vector<NodeType::Pointer> ConditionNodes2 (3);
//             
//             ConditionNodes2[0] = pNode1;
//             ConditionNodes2[1] = pNode3;
//             ConditionNodes2[2] = pNode4;
//             
//             Triangle2D3 <Node<3>> triangle2( ConditionNodes2 );
//             
//             // We calculate the integral of the mass matrix (assuming constant density)
//             GeometryNodeType::IntegrationPointsArrayType IntegrationPointsQuadrilateral = Quadrature<QuadrilateralGaussLegendreIntegrationPoints2, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
//             GeometryNodeType::IntegrationPointsArrayType IntegrationPointsTriangle = Quadrature<TriangleGaussLegendreIntegrationPoints5, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
//             
//             boost::numeric::ublas::bounded_matrix<double, 4, 4> MassMatrix0 = ZeroMatrix(4, 4);
//             
//             for (unsigned int PointNumber = 0; PointNumber < IntegrationPointsQuadrilateral.size(); PointNumber++)
//             {
//                 Vector N;
//                 const PointType& LocalPoint = IntegrationPointsQuadrilateral[PointNumber].Coordinates();
//                 quadrilateral0.ShapeFunctionsValues( N, LocalPoint );
//                 const double DetJ = quadrilateral0.DeterminantOfJacobian( LocalPoint );
//                 const double weight = IntegrationPointsQuadrilateral[PointNumber].Weight();
//                 
//                 for (unsigned int inode = 0; inode < 4; inode++)
//                 {
//                     for (unsigned int jnode = 0; jnode < 4; jnode++)
//                     {
//                         MassMatrix0(inode, jnode) += DetJ * weight * N[inode] * N[jnode];
//                     }
//                 }
//             }
//      
//             boost::numeric::ublas::bounded_matrix<double, 4, 4> MassMatrix1 = ZeroMatrix(4, 4);
//             
//             array_1d<Point<3>, 3> LocalPointsList1;
//             array_1d<Point<3>, 3> LocalPointsList2;
//             
//             for (unsigned int i_point = 0; i_point < 3; i_point++)
//             {
//                 PointType LocalPoint;
//                 quadrilateral0.PointLocalCoordinates(LocalPoint, triangle1[i_point]);
//                 LocalPointsList1[i_point] = LocalPoint;
//                 quadrilateral0.PointLocalCoordinates(LocalPoint, triangle2[i_point]);
//                 LocalPointsList2[i_point] = LocalPoint;
//             }
//             
//             for (unsigned int PointNumber = 0; PointNumber < IntegrationPointsTriangle.size(); PointNumber++)
//             {
//                 Vector Naux1;
//                 Vector Naux2;
//                 
//                 Vector N1;
//                 Vector N2;
//                 
//                 const PointType& LocalPoint0 = IntegrationPointsTriangle[PointNumber].Coordinates();
//             
//                 triangle1.ShapeFunctionsValues( Naux1, LocalPoint0 );
//                 triangle2.ShapeFunctionsValues( Naux2, LocalPoint0 );
//                 
//                 PointType LocalPoint1;
//                 PointType LocalPoint2; 
//                 
//                 LocalPoint1.Coordinates() = ZeroVector(3);
//                 LocalPoint2.Coordinates() = ZeroVector(3);
//                 
//                 for (unsigned iCoor = 1; iCoor <= 2; iCoor++)
//                 {
//                     for (unsigned int iN = 0; iN < 3; iN++)
//                     {
//                         LocalPoint1.Coordinate(iCoor) += Naux1[iN] * LocalPointsList1[iN].Coordinate(iCoor);
//                         LocalPoint2.Coordinate(iCoor) += Naux2[iN] * LocalPointsList2[iN].Coordinate(iCoor);
//                     }
//                 }
//                 
//                 PointType gp_global;
//                 
//                 triangle1.GlobalCoordinates(gp_global, LocalPoint0);
//                 PointType LocalPoint1aux;
//                 LocalPoint1aux.Coordinates() = ZeroVector(3);
//                 quadrilateral0.PointLocalCoordinates(LocalPoint1aux, gp_global);
//                 
//                 triangle2.GlobalCoordinates(gp_global, LocalPoint0);
//                 PointType LocalPoint2aux; 
//                 LocalPoint2aux.Coordinates() = ZeroVector(3);
//                 quadrilateral0.PointLocalCoordinates(LocalPoint2aux, gp_global);
//                 
// //                 std::cout << std::endl;
// //                 KRATOS_WATCH(LocalPoint1);
// //                 KRATOS_WATCH(LocalPoint1aux);
// //                 KRATOS_WATCH(LocalPoint2);
// //                 KRATOS_WATCH(LocalPoint2aux);
//                 
//                 quadrilateral0.ShapeFunctionsValues( N1, LocalPoint1 );
//                 quadrilateral0.ShapeFunctionsValues( N2, LocalPoint2 );
//                 
//                 const double DetJ1  = triangle1.DeterminantOfJacobian( LocalPoint0 );
//                 const double DetJ2  = triangle2.DeterminantOfJacobian( LocalPoint0 );
//                 
//                 const double weight = IntegrationPointsTriangle[PointNumber].Weight();
//                 
//                 for (unsigned int inode = 0; inode < 4; inode++)
//                 {
//                     for (unsigned int jnode = 0; jnode < 4; jnode++)
//                     {                        
//                         MassMatrix1(inode, jnode ) += DetJ1 * weight * N1[inode] * N1[jnode] 
//                                                    +  DetJ2 * weight * N2[inode] * N2[jnode];
//                     }
//                 }
//             }
// 
//             
//             // Debug
// //             KRATOS_WATCH(MassMatrix0)
// //             KRATOS_WATCH(MassMatrix1)
//             
//             const double Tolerance = 1.0e-6;
//             for (unsigned int inode = 0; inode < 4; inode++)
//             {
//                 for (unsigned int jnode = 0; jnode < 4; jnode++)
//                 {
//                     KRATOS_CHECK_NEAR(MassMatrix0(inode,jnode), MassMatrix1(inode,jnode), Tolerance);
//                 }
//             }
//         }
        
    } // namespace Testing
}  // namespace Kratos.
