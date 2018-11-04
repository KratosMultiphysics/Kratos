// KRATOS  ___|  |       |       |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//           | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License: BSD License
//   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Vicente Mataix Ferrandiz
//

// System includes
#include <set>

// External includes

// Project includes
#include "testing/testing.h"
#include "spaces/ublas_space.h"
#include "includes/properties.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "utilities/math_utils.h"

/* GAUSS-LEGENDRE */
#include "integration/triangle_gauss_legendre_integration_points.h"
#include "integration/quadrilateral_gauss_legendre_integration_points.h"

/* Utilities */
#include "utilities/mortar_utilities.h"
#include "utilities/exact_mortar_segmentation_utility.h"

namespace Kratos 
{
    namespace Testing 
    {

        typedef Point                                                     PointType;
        typedef Node<3>                                                    NodeType;
        typedef Geometry<NodeType>                                 GeometryNodeType;
        typedef Geometry<PointType>                               GeometryPointType;
        typedef std::size_t                                               IndexType;
        
        ///Type definition for integration methods
        typedef GeometryData::IntegrationMethod                   IntegrationMethod;
        typedef IntegrationPoint<2>                            IntegrationPointType;
        typedef GeometryNodeType::IntegrationPointsArrayType integration_pointsType;

        /** 
         * Checks if the criteria for computing the integral is the correct one
         * Checks mass matrix computed
         */
    
        KRATOS_TEST_CASE_IN_SUITE(MassMatrixIntegrationTriangle, KratosContactStructuralMechanicsFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
            
            // First we create the nodes 
            NodeType::Pointer p_node_1 = r_model_part.CreateNewNode(0,-0.2,0.1,0.0);
            NodeType::Pointer p_node_2 = r_model_part.CreateNewNode(1,1.0,0.1,0.0);
            NodeType::Pointer p_node_3 = r_model_part.CreateNewNode(2,0.2,1.2,0.0);
            NodeType::Pointer p_node_4 = r_model_part.CreateNewNode(3,0.6,0.4,0.0);
            
            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (3);
            
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            condition_nodes_0[2] = p_node_3;
            
            Triangle3D3 <NodeType> triangle0( PointerVector<NodeType>{condition_nodes_0} );
            
            std::vector<NodeType::Pointer> condition_nodes_1 (3);
            
            condition_nodes_1[0] = p_node_1;
            condition_nodes_1[1] = p_node_2;
            condition_nodes_1[2] = p_node_4;
            
            Triangle3D3 <NodeType> triangle_1( PointerVector<NodeType>{condition_nodes_1} );
            
            std::vector<NodeType::Pointer> condition_nodes_2 (3);
            
            condition_nodes_2[0] = p_node_2;
            condition_nodes_2[1] = p_node_3;
            condition_nodes_2[2] = p_node_4;
            
            Triangle3D3 <NodeType> triangle_2( PointerVector<NodeType>{condition_nodes_2} );
            
            std::vector<NodeType::Pointer> condition_nodes_3 (3);
            
            condition_nodes_3[0] = p_node_3;
            condition_nodes_3[1] = p_node_1;
            condition_nodes_3[2] = p_node_4;
            
            Triangle3D3 <NodeType> triangle_3( PointerVector<NodeType>{condition_nodes_3} );
            
            // We calculate the integral of the mass matrix (assuming constant density)
            GeometryNodeType::IntegrationPointsArrayType integration_points = Quadrature<TriangleGaussLegendreIntegrationPoints2, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            
            BoundedMatrix<double, 3, 3> mass_matrix_0 = ZeroMatrix(3, 3);
            
            for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number)
            {
                Vector N;
                const PointType& local_point = integration_points[point_number].Coordinates();
                triangle0.ShapeFunctionsValues( N, local_point );
                const double det_j = triangle0.DeterminantOfJacobian( local_point );
                const double weight = integration_points[point_number].Weight();
                
                for (IndexType i_node = 0; i_node < 3; ++i_node)
                    for (IndexType j_node = 0; j_node < 3; ++j_node)
                        mass_matrix_0(i_node, j_node) += det_j * weight * N(i_node) * N(j_node);
            }
     
            BoundedMatrix<double, 3, 3> mass_matrix_1 = ZeroMatrix(3, 3);
            
            for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number)
            {
                Vector N1, N2, N3;
                
                const PointType& local_point_0 = integration_points[point_number].Coordinates();
                
                PointType gp_global;
                
                triangle_1.GlobalCoordinates(gp_global, local_point_0);
                PointType local_point_1;
                triangle0.PointLocalCoordinates(local_point_1, gp_global);
                
                triangle_2.GlobalCoordinates(gp_global, local_point_0);
                PointType local_point_2; 
                triangle0.PointLocalCoordinates(local_point_2, gp_global);
                
                triangle_3.GlobalCoordinates(gp_global, local_point_0);
                PointType local_point_3; 
                triangle0.PointLocalCoordinates(local_point_3, gp_global);
                
                triangle0.ShapeFunctionsValues( N1, local_point_1 );
                triangle0.ShapeFunctionsValues( N2, local_point_2 );
                triangle0.ShapeFunctionsValues( N3, local_point_3 );
                
                const double det_j_1 = triangle_1.DeterminantOfJacobian( local_point_0 );
                const double det_j_2 = triangle_2.DeterminantOfJacobian( local_point_0 );
                const double det_j_3 = triangle_3.DeterminantOfJacobian( local_point_0 );
                
                const double weight = integration_points[point_number].Weight();
                
                for (IndexType i_node = 0; i_node < 3; ++i_node)
                {
                    for (IndexType j_node = 0; j_node < 3; ++j_node)
                    {
                        mass_matrix_1(i_node, j_node) += det_j_1 * weight * N1[i_node] * N1[j_node] \
                                                       + det_j_2 * weight * N2[i_node] * N2[j_node] \
                                                       + det_j_3 * weight * N3[i_node] * N3[j_node];
                    }
                }
            }
            
            const double tolerance = 1.0e-6;
            for (IndexType i_node = 0; i_node < 3; ++i_node)
                for (IndexType j_node = 0; j_node < 3; ++j_node)
                    KRATOS_CHECK_NEAR(mass_matrix_0(i_node,j_node), mass_matrix_1(i_node,j_node), tolerance);
        }
        
        /** 
         * Checks if the criteria for computing the integral is the correct one
         * Checks mass matrix computed
         */
    
        KRATOS_TEST_CASE_IN_SUITE(MassMatrixIntegrationQuadrilateral, KratosContactStructuralMechanicsFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
            
            // First we create the nodes 
            NodeType::Pointer p_node_1 = r_model_part.CreateNewNode(0,   0.0,  0.0, 0.0);
            NodeType::Pointer p_node_2 = r_model_part.CreateNewNode(1,   1.0,- 0.1, 0.0);
            NodeType::Pointer p_node_3 = r_model_part.CreateNewNode(2,   1.2,  1.1, 0.0);
            NodeType::Pointer p_node_4 = r_model_part.CreateNewNode(3, - 0.1,  1.3, 0.0);
            
            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (4);
            
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            condition_nodes_0[2] = p_node_3;
            condition_nodes_0[3] = p_node_4;
            
            Quadrilateral3D4 <NodeType> quadrilateral_0( PointerVector<NodeType>{condition_nodes_0} );
            
            std::vector<NodeType::Pointer> condition_nodes_1 (3);
            
            condition_nodes_1[0] = p_node_1;
            condition_nodes_1[1] = p_node_2;
            condition_nodes_1[2] = p_node_3;
            
            Triangle3D3 <NodeType> triangle_1( PointerVector<NodeType>{condition_nodes_1} );
            
            std::vector<NodeType::Pointer> condition_nodes_2 (3);
            
            condition_nodes_2[0] = p_node_1;
            condition_nodes_2[1] = p_node_3;
            condition_nodes_2[2] = p_node_4;
            
            Triangle3D3 <NodeType> triangle_2( PointerVector<NodeType>{condition_nodes_2} );
            
            // We calculate the integral of the mass matrix (assuming constant density)
            GeometryNodeType::IntegrationPointsArrayType integration_pointsQuadrilateral = Quadrature<QuadrilateralGaussLegendreIntegrationPoints2, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            GeometryNodeType::IntegrationPointsArrayType integration_pointsTriangle = Quadrature<TriangleGaussLegendreIntegrationPoints5, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            
            BoundedMatrix<double, 4, 4> mass_matrix_0 = ZeroMatrix(4, 4);
            
            for (IndexType point_number = 0; point_number < integration_pointsQuadrilateral.size(); ++point_number)
            {
                Vector N;
                const PointType& local_point = integration_pointsQuadrilateral[point_number].Coordinates();
                quadrilateral_0.ShapeFunctionsValues( N, local_point );
                const double det_j = quadrilateral_0.DeterminantOfJacobian( local_point );
                const double weight = integration_pointsQuadrilateral[point_number].Weight();
                
                for (IndexType i_node = 0; i_node < 4; ++i_node)
                    for (IndexType j_node = 0; j_node < 4; ++j_node)
                        mass_matrix_0(i_node, j_node) += det_j * weight * N[i_node] * N[j_node];
            }
     
            BoundedMatrix<double, 4, 4> mass_matrix_1 = ZeroMatrix(4, 4);
            
            for (IndexType point_number = 0; point_number < integration_pointsTriangle.size(); ++point_number)
            {
                Vector N1, N2;
                
                const PointType& local_point_0 = integration_pointsTriangle[point_number].Coordinates();
                
                PointType gp_global;
                
                triangle_1.GlobalCoordinates(gp_global, local_point_0);
                PointType local_point_1;
                quadrilateral_0.PointLocalCoordinates(local_point_1, gp_global);
                
                triangle_2.GlobalCoordinates(gp_global, local_point_0);
                PointType local_point_2; 
                quadrilateral_0.PointLocalCoordinates(local_point_2, gp_global);
                
                quadrilateral_0.ShapeFunctionsValues( N1, local_point_1 );
                quadrilateral_0.ShapeFunctionsValues( N2, local_point_2 );
                
                const double det_j_1  = triangle_1.DeterminantOfJacobian( local_point_0 );
                const double det_j_2  = triangle_2.DeterminantOfJacobian( local_point_0 );
                
                const double weight = integration_pointsTriangle[point_number].Weight();
                
                for (IndexType i_node = 0; i_node < 4; ++i_node)
                {
                    for (IndexType j_node = 0; j_node < 4; ++j_node)
                    {                        
                        mass_matrix_1(i_node, j_node ) += det_j_1 * weight * N1[i_node] * N1[j_node] 
                                                        + det_j_2 * weight * N2[i_node] * N2[j_node];
                    }
                }
            }

            
            // Debug
//             KRATOS_WATCH(mass_matrix_0)
//             KRATOS_WATCH(mass_matrix_1)
            
            const double tolerance = 1.0e-6;
            for (IndexType i_node = 0; i_node < 4; ++i_node)
                for (IndexType j_node = 0; j_node < 4; ++j_node)
                    KRATOS_CHECK_NEAR(mass_matrix_0(i_node,j_node), mass_matrix_1(i_node,j_node), tolerance);
        }
        
        /** 
         * Checks if the criteria for computing the integral is the correct one
         * Checks mass matrix computed
         */
    
        KRATOS_TEST_CASE_IN_SUITE(MassMatrixIntegrationQuadrilateralDeformed, KratosContactStructuralMechanicsFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
            
            // First we create the nodes 
            NodeType::Pointer p_node_0 = r_model_part.CreateNewNode(0,   0.5,  0.4, 0.0);
            NodeType::Pointer p_node_1 = r_model_part.CreateNewNode(1,   0.0,  0.0, 0.0);
            NodeType::Pointer p_node_2 = r_model_part.CreateNewNode(2,   1.0,- 0.1, 0.0);
            NodeType::Pointer p_node_3 = r_model_part.CreateNewNode(3,   1.2,  1.1, 0.0);
            NodeType::Pointer p_node_4 = r_model_part.CreateNewNode(4, - 0.1,  1.3, 0.0);
            
            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (4);
            
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            condition_nodes_0[2] = p_node_3;
            condition_nodes_0[3] = p_node_4;
            
            Quadrilateral3D4 <NodeType> quadrilateral_0( PointerVector<NodeType>{condition_nodes_0} );
            
            std::vector<NodeType::Pointer> condition_nodes_1 (3);
            
            condition_nodes_1[0] = p_node_1;
            condition_nodes_1[1] = p_node_2;
            condition_nodes_1[2] = p_node_0;
            
            Triangle3D3 <NodeType> triangle_1( PointerVector<NodeType>{condition_nodes_1} );
            
            std::vector<NodeType::Pointer> condition_nodes_2 (3);
            
            condition_nodes_2[0] = p_node_2;
            condition_nodes_2[1] = p_node_3;
            condition_nodes_2[2] = p_node_0;
            
            Triangle3D3 <NodeType> triangle_2( PointerVector<NodeType>{condition_nodes_2} );
            
            std::vector<NodeType::Pointer> condition_nodes_3 (3);
            
            condition_nodes_3[0] = p_node_3;
            condition_nodes_3[1] = p_node_4;
            condition_nodes_3[2] = p_node_0;
            
            Triangle3D3 <NodeType> triangle_3( PointerVector<NodeType>{condition_nodes_3} );
            
            std::vector<NodeType::Pointer> condition_nodes_4 (3);
            
            condition_nodes_4[0] = p_node_4;
            condition_nodes_4[1] = p_node_1;
            condition_nodes_4[2] = p_node_0;
            
            Triangle3D3 <NodeType> triangle_4( PointerVector<NodeType>{condition_nodes_4} );
            
            // We calculate the integral of the mass matrix (assuming constant density)
            GeometryNodeType::IntegrationPointsArrayType integration_pointsQuadrilateral = Quadrature<QuadrilateralGaussLegendreIntegrationPoints2, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            GeometryNodeType::IntegrationPointsArrayType integration_pointsTriangle = Quadrature<TriangleGaussLegendreIntegrationPoints5, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            
            BoundedMatrix<double, 4, 4> mass_matrix_0 = ZeroMatrix(4, 4);
            
            for (IndexType point_number = 0; point_number < integration_pointsQuadrilateral.size(); ++point_number)
            {
                Vector N;
                const PointType& local_point = integration_pointsQuadrilateral[point_number].Coordinates();
                quadrilateral_0.ShapeFunctionsValues( N, local_point );
                const double det_j = quadrilateral_0.DeterminantOfJacobian( local_point );
                const double weight = integration_pointsQuadrilateral[point_number].Weight();
                
                for (IndexType i_node = 0; i_node < 4; ++i_node)
                    for (IndexType j_node = 0; j_node < 4; ++j_node)
                        mass_matrix_0(i_node, j_node) += det_j * weight * N[i_node] * N[j_node];
            }
     
            BoundedMatrix<double, 4, 4> mass_matrix_1 = ZeroMatrix(4, 4);
            
            for (IndexType point_number = 0; point_number < integration_pointsTriangle.size(); ++point_number)
            {
                Vector N1, N2, N3, N4;
                
                const PointType& local_point_0 = integration_pointsTriangle[point_number].Coordinates();
                
                PointType gp_global;
                
                triangle_1.GlobalCoordinates(gp_global, local_point_0);
                PointType local_point_1;
                quadrilateral_0.PointLocalCoordinates(local_point_1, gp_global);
                
                triangle_2.GlobalCoordinates(gp_global, local_point_0);
                PointType local_point_2; 
                quadrilateral_0.PointLocalCoordinates(local_point_2, gp_global);
                
                triangle_3.GlobalCoordinates(gp_global, local_point_0);
                PointType local_point_3; 
                quadrilateral_0.PointLocalCoordinates(local_point_3, gp_global);
                
                triangle_4.GlobalCoordinates(gp_global, local_point_0);
                PointType local_point_4; 
                quadrilateral_0.PointLocalCoordinates(local_point_4, gp_global);
                
                quadrilateral_0.ShapeFunctionsValues( N1, local_point_1 );
                quadrilateral_0.ShapeFunctionsValues( N2, local_point_2 );
                quadrilateral_0.ShapeFunctionsValues( N3, local_point_3 );
                quadrilateral_0.ShapeFunctionsValues( N4, local_point_4 );
                
                const double det_j_1  = triangle_1.DeterminantOfJacobian( local_point_0 );
                const double det_j_2  = triangle_2.DeterminantOfJacobian( local_point_0 );
                const double det_j_3  = triangle_3.DeterminantOfJacobian( local_point_0 );
                const double det_j_4  = triangle_4.DeterminantOfJacobian( local_point_0 );
                
                const double weight = integration_pointsTriangle[point_number].Weight();
                
                for (IndexType i_node = 0; i_node < 4; ++i_node)
                {
                    for (IndexType j_node = 0; j_node < 4; ++j_node)
                    {                        
                        mass_matrix_1(i_node, j_node ) += det_j_1 * weight * N1[i_node] * N1[j_node] 
                                                       +  det_j_2 * weight * N2[i_node] * N2[j_node]
                                                       +  det_j_3 * weight * N3[i_node] * N3[j_node]
                                                       +  det_j_4 * weight * N4[i_node] * N4[j_node];
                    }
                }
            }
                        
//             // Debug
//             KRATOS_WATCH(mass_matrix_0)
//             KRATOS_WATCH(mass_matrix_1)
            
            const double tolerance = 1.0e-6;
            for (IndexType i_node = 0; i_node < 4; ++i_node)
                for (IndexType j_node = 0; j_node < 4; ++j_node)
                    KRATOS_CHECK_NEAR(mass_matrix_0(i_node,j_node), mass_matrix_1(i_node,j_node), tolerance);
            
            array_1d<double, 3> disp_array = ZeroVector(3);
            disp_array[0] = 0.1;
            disp_array[1] = 0.2;
            
            p_node_0->Coordinates() += 0.5  * disp_array;
            p_node_1->Coordinates() += 1.0  * disp_array;
            p_node_2->Coordinates() += 1.5  * disp_array;
            p_node_3->Coordinates() += 0.75 * disp_array;
            p_node_4->Coordinates() += 1.0  * disp_array;
            
            mass_matrix_0 = ZeroMatrix(4, 4);
            
            for (IndexType point_number = 0; point_number < integration_pointsQuadrilateral.size(); ++point_number)
            {
                Vector N;
                const PointType& local_point = integration_pointsQuadrilateral[point_number].Coordinates();
                quadrilateral_0.ShapeFunctionsValues( N, local_point );
                const double det_j = quadrilateral_0.DeterminantOfJacobian( local_point );
                const double weight = integration_pointsQuadrilateral[point_number].Weight();
                
                for (IndexType i_node = 0; i_node < 4; ++i_node)
                    for (IndexType j_node = 0; j_node < 4; ++j_node)
                        mass_matrix_0(i_node, j_node) += det_j * weight * N[i_node] * N[j_node];
            }
            
            mass_matrix_1 = ZeroMatrix(4, 4);
            
            for (IndexType point_number = 0; point_number < integration_pointsTriangle.size(); ++point_number)
            {
                Vector N1, N2, N3, N4;
                
                const PointType& local_point_0 = integration_pointsTriangle[point_number].Coordinates();
                
                PointType gp_global;
                
                triangle_1.GlobalCoordinates(gp_global, local_point_0);
                PointType local_point_1;
                quadrilateral_0.PointLocalCoordinates(local_point_1, gp_global);
                
                triangle_2.GlobalCoordinates(gp_global, local_point_0);
                PointType local_point_2; 
                quadrilateral_0.PointLocalCoordinates(local_point_2, gp_global);
                
                triangle_3.GlobalCoordinates(gp_global, local_point_0);
                PointType local_point_3; 
                quadrilateral_0.PointLocalCoordinates(local_point_3, gp_global);
                
                triangle_4.GlobalCoordinates(gp_global, local_point_0);
                PointType local_point_4; 
                quadrilateral_0.PointLocalCoordinates(local_point_4, gp_global);
                
                quadrilateral_0.ShapeFunctionsValues( N1, local_point_1 );
                quadrilateral_0.ShapeFunctionsValues( N2, local_point_2 );
                quadrilateral_0.ShapeFunctionsValues( N3, local_point_3 );
                quadrilateral_0.ShapeFunctionsValues( N4, local_point_4 );
                
                const double det_j_1 = triangle_1.DeterminantOfJacobian( local_point_0 );
                const double det_j_2 = triangle_2.DeterminantOfJacobian( local_point_0 );
                const double det_j_3 = triangle_3.DeterminantOfJacobian( local_point_0 );
                const double det_j_4 = triangle_4.DeterminantOfJacobian( local_point_0 );
                
                const double weight = integration_pointsTriangle[point_number].Weight();
                
                for (IndexType i_node = 0; i_node < 4; ++i_node)
                {
                    for (IndexType j_node = 0; j_node < 4; ++j_node)
                    {                        
                        mass_matrix_1(i_node, j_node ) += det_j_1 * weight * N1[i_node] * N1[j_node] 
                                                       +  det_j_2 * weight * N2[i_node] * N2[j_node]
                                                       +  det_j_3 * weight * N3[i_node] * N3[j_node]
                                                       +  det_j_4 * weight * N4[i_node] * N4[j_node];
                    }
                }
            }
                        
//             // Debug
//             KRATOS_WATCH(mass_matrix_0)
//             KRATOS_WATCH(mass_matrix_1)
            
            for (IndexType i_node = 0; i_node < 4; ++i_node)
                for (IndexType j_node = 0; j_node < 4; ++j_node)
                    KRATOS_CHECK_NEAR(mass_matrix_0(i_node,j_node), mass_matrix_1(i_node,j_node), tolerance);
            
        }
        
        /** 
         * Checks if the criteria for computing the integral is the correct one interpolating the local coordinates
         * Checks mass matrix computed
         */
        
        KRATOS_TEST_CASE_IN_SUITE(TestCheckRotation, KratosContactStructuralMechanicsFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
            
            // First we create the nodes 
            NodeType::Pointer p_node_1 = r_model_part.CreateNewNode(0,   0.0,  0.0, 0.1);
            NodeType::Pointer p_node_2 = r_model_part.CreateNewNode(1,   1.0,- 0.1, 0.0);
            NodeType::Pointer p_node_3 = r_model_part.CreateNewNode(2,   1.2,  1.1, 0.2);
            
            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (3);
            
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            condition_nodes_0[2] = p_node_3;
            
            Triangle3D3 <NodeType> triangle_0( PointerVector<NodeType>{condition_nodes_0} );
            
            // We define the condition tangents
            const array_1d<double, 3> slave_tangent_xi  = (triangle_0[1].Coordinates() - triangle_0[0].Coordinates())/norm_2(triangle_0[1].Coordinates() - triangle_0[0].Coordinates());
            array_1d<double, 3> aux_coords;
            triangle_0.PointLocalCoordinates(aux_coords, triangle_0.Center());
            const array_1d<double, 3> normal = triangle_0.UnitNormal(aux_coords);
            array_1d<double, 3> slave_tangent_eta;
            MathUtils<double>::CrossProduct(slave_tangent_eta, normal, slave_tangent_xi);
            
            // We define the auxiliar geometry
            std::vector<PointType::Pointer> points_array  (3);
            for (IndexType i_node = 0; i_node < 3; ++i_node)
            {
                PointType aux_point;
                aux_point.Coordinates() = triangle_0[i_node].Coordinates();
                points_array[i_node] = PointType::Pointer( new PointType(aux_point) );
            }
            
            Triangle3D3 <PointType> aux_geometry(  PointerVector<PointType>{points_array}  );
            const PointType center = aux_geometry.Center();
            
            // Before clipping we rotate to a XY plane
            for (IndexType i_node = 0; i_node < 3; ++i_node)
            {
                MortarUtilities::RotatePoint( aux_geometry[i_node], center, slave_tangent_xi, slave_tangent_eta, false);
                MortarUtilities::RotatePoint( aux_geometry[i_node], center, slave_tangent_xi, slave_tangent_eta, true);
            }
            
            const double tolerance = 1.0e-6;
            for (IndexType i_node = 0; i_node < 3; ++i_node)
            {
                const array_1d<double, 3>& coords1 = aux_geometry[i_node].Coordinates();
                const array_1d<double, 3>& coords2 = triangle_0[i_node].Coordinates();
                for (IndexType jdim = 0; jdim < 3; jdim++)
                    KRATOS_CHECK_NEAR(coords1[jdim], coords2[jdim], tolerance);
            }
        }
        
    } // namespace Testing
}  // namespace Kratos.
