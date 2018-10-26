//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:   BSD License
//      Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "testing/testing.h"
#include "includes/model_part.h"

/* Utilities */
#include "utilities/binbased_fast_point_locator.h"
#include "utilities/binbased_fast_point_locator_conditions.h"

namespace Kratos 
{
    namespace Testing 
    {
        typedef Node<3> NodeType;
        
        /** 
        * Checks the correct work of the binbased fast point locator
        * Test triangle 
        */

        KRATOS_TEST_CASE_IN_SUITE(TestBinBasedFastPointLocator1, KratosCoreFastSuite)
        {
            Model current_model;
            ModelPart& this_model_part = current_model.CreateModelPart("Main");
            this_model_part.SetBufferSize(2);
            
            Properties::Pointer p_elem_prop = this_model_part.pGetProperties(0);
            
            auto& process_info = this_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;
            
            // First we create the nodes 
            NodeType::Pointer p_node_1 = this_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_2 = this_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_3 = this_model_part.CreateNewNode(3, 1.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_4 = this_model_part.CreateNewNode(4, 0.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_5 = this_model_part.CreateNewNode(5, 2.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_6 = this_model_part.CreateNewNode(6, 2.0 , 1.0 , 0.0);
            
            // Now we create the "conditions"
            std::vector<NodeType::Pointer> element_nodes_0 (3);
            element_nodes_0[0] = p_node_1;
            element_nodes_0[1] = p_node_2;
            element_nodes_0[2] = p_node_3;
            Triangle2D3 <NodeType> triangle_0( PointerVector<NodeType>{element_nodes_0} );
            
            std::vector<NodeType::Pointer> element_nodes_1 (3);
            element_nodes_1[0] = p_node_1;
            element_nodes_1[1] = p_node_3;
            element_nodes_1[2] = p_node_4;
            Triangle2D3 <NodeType> triangle_1( PointerVector<NodeType>{element_nodes_1} );
            
            std::vector<NodeType::Pointer> element_nodes_2 (3);
            element_nodes_2[0] = p_node_2;
            element_nodes_2[1] = p_node_5;
            element_nodes_2[2] = p_node_3;
            Triangle2D3 <NodeType> triangle_2( PointerVector<NodeType>{element_nodes_2} );
            
            std::vector<NodeType::Pointer> element_nodes_3 (3);
            element_nodes_3[0] = p_node_5;
            element_nodes_3[1] = p_node_6;
            element_nodes_3[2] = p_node_3;
            Triangle2D3 <NodeType> triangle_3( PointerVector<NodeType>{element_nodes_3} );
            
            Element::Pointer p_elem_0 = this_model_part.CreateNewElement("Element2D3N", 1, triangle_0, p_elem_prop);
            Element::Pointer p_elem_1 = this_model_part.CreateNewElement("Element2D3N", 2, triangle_1, p_elem_prop);
            Element::Pointer p_elem_2 = this_model_part.CreateNewElement("Element2D3N", 3, triangle_2, p_elem_prop);
            Element::Pointer p_elem_3 = this_model_part.CreateNewElement("Element2D3N", 4, triangle_3, p_elem_prop);
            
            // We create the locator
            auto point_locator = BinBasedFastPointLocator<2>(this_model_part);
            point_locator.UpdateSearchDatabase();

            array_1d<double, 3> coordinates(3, 0.0);
            coordinates[0] = 0.5;
            coordinates[1] = 0.5;
            Vector shape_functions;
            Element::Pointer p_element;
            bool is_found = point_locator.FindPointOnMeshSimplified(coordinates, shape_functions, p_element, 1000, 5.0e-2);

            Vector ref_shape_functions(3);
            ref_shape_functions[0] = 0.5;
            ref_shape_functions[1] = 0.0;
            ref_shape_functions[2] = 0.5;

            const double tolerance = 1.0e-16;
            KRATOS_CHECK(is_found);
            KRATOS_CHECK_EQUAL(p_element->Id(), p_elem_0->Id());
            KRATOS_CHECK_LESS_EQUAL(norm_2(shape_functions - ref_shape_functions), tolerance);

            coordinates[0] = -0.5;
            coordinates[1] = -0.5;
            is_found = point_locator.FindPointOnMeshSimplified(coordinates, shape_functions, p_element, 1000, 5.0e-2);
            KRATOS_CHECK_IS_FALSE(is_found);
        }
        
        /** 
        * Checks the correct work of the binbased fast point locator
        * Test tetrahedra
        */

        KRATOS_TEST_CASE_IN_SUITE(TestBinBasedFastPointLocator2, KratosCoreFastSuite)
        {
            Model current_model;
            ModelPart& this_model_part = current_model.CreateModelPart("Main");
            this_model_part.SetBufferSize(2);
            
            Properties::Pointer p_elem_prop = this_model_part.pGetProperties(0);
            
            auto& process_info = this_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;
            
            // First we create the nodes 
            NodeType::Pointer p_node_1 = this_model_part.CreateNewNode(1 , 0.0 , 1.0 , 1.0);
            NodeType::Pointer p_node_2 = this_model_part.CreateNewNode(2 , 0.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_3 = this_model_part.CreateNewNode(3 , 0.0 , 0.0 , 1.0);
            NodeType::Pointer p_node_4 = this_model_part.CreateNewNode(4 , 1.0 , 1.0 , 1.0);
            NodeType::Pointer p_node_5 = this_model_part.CreateNewNode(5 , 0.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_6 = this_model_part.CreateNewNode(6 , 1.0 , 1.0 , 0.0);
            
            NodeType::Pointer p_node_7 = this_model_part.CreateNewNode(7 , 1.0 , 0.0 , 1.0);
            NodeType::Pointer p_node_8 = this_model_part.CreateNewNode(8 , 1.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_9 = this_model_part.CreateNewNode(9 , 2.0 , 1.0 , 1.0);
            NodeType::Pointer p_node_10 = this_model_part.CreateNewNode(10 , 2.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_11 = this_model_part.CreateNewNode(11 , 2.0 , 0.0 , 1.0);
            NodeType::Pointer p_node_12 = this_model_part.CreateNewNode(12 , 2.0 , 0.0 , 0.0);
            
            // Now we create the "conditions"
            std::vector<NodeType::Pointer> element_nodes_0 (4);
            element_nodes_0[0] = p_node_12;
            element_nodes_0[1] = p_node_10;
            element_nodes_0[2] = p_node_8;
            element_nodes_0[3] = p_node_9;
            Tetrahedra3D4 <NodeType> tetrahedra_0( PointerVector<NodeType>{element_nodes_0} );
            
            std::vector<NodeType::Pointer> element_nodes_1 (4);
            element_nodes_1[0] = p_node_4;
            element_nodes_1[1] = p_node_6;
            element_nodes_1[2] = p_node_9;
            element_nodes_1[3] = p_node_7;
            Tetrahedra3D4 <NodeType> tetrahedra_1( PointerVector<NodeType>{element_nodes_1} );
            
            std::vector<NodeType::Pointer> element_nodes_2 (4);
            element_nodes_2[0] = p_node_11;
            element_nodes_2[1] = p_node_7;
            element_nodes_2[2] = p_node_9;
            element_nodes_2[3] = p_node_8;
            Tetrahedra3D4 <NodeType> tetrahedra_2( PointerVector<NodeType>{element_nodes_2} );
            
            std::vector<NodeType::Pointer> element_nodes_3 (4);
            element_nodes_3[0] = p_node_5;
            element_nodes_3[1] = p_node_3;
            element_nodes_3[2] = p_node_8;
            element_nodes_3[3] = p_node_6;
            Tetrahedra3D4 <NodeType> tetrahedra_3( PointerVector<NodeType>{element_nodes_3} );
            
            std::vector<NodeType::Pointer> element_nodes_4 (4);
            element_nodes_4[0] = p_node_4;
            element_nodes_4[1] = p_node_6;
            element_nodes_4[2] = p_node_7;
            element_nodes_4[3] = p_node_3;
            Tetrahedra3D4 <NodeType> tetrahedra_4( PointerVector<NodeType>{element_nodes_4} );
            
            std::vector<NodeType::Pointer> element_nodes_5 (4);
            element_nodes_5[0] = p_node_2;
            element_nodes_5[1] = p_node_3;
            element_nodes_5[2] = p_node_5;
            element_nodes_5[3] = p_node_6;
            Tetrahedra3D4 <NodeType> tetrahedra_5( PointerVector<NodeType>{element_nodes_5} );
            
            std::vector<NodeType::Pointer> element_nodes_6 (4);
            element_nodes_6[0] = p_node_10;
            element_nodes_6[1] = p_node_9;
            element_nodes_6[2] = p_node_6;
            element_nodes_6[3] = p_node_8;
            Tetrahedra3D4 <NodeType> tetrahedra_6( PointerVector<NodeType>{element_nodes_6} );
            
            std::vector<NodeType::Pointer> element_nodes_7 (4);
            element_nodes_7[0] = p_node_7;
            element_nodes_7[1] = p_node_8;
            element_nodes_7[2] = p_node_3;
            element_nodes_7[3] = p_node_6;
            Tetrahedra3D4 <NodeType> tetrahedra_7( PointerVector<NodeType>{element_nodes_7} );
            
            std::vector<NodeType::Pointer> element_nodes_8 (4);
            element_nodes_8[0] = p_node_7;
            element_nodes_8[1] = p_node_8;
            element_nodes_8[2] = p_node_6;
            element_nodes_8[3] = p_node_9;
            Tetrahedra3D4 <NodeType> tetrahedra_8( PointerVector<NodeType>{element_nodes_8} );
            
            std::vector<NodeType::Pointer> element_nodes_9 (4);
            element_nodes_9[0] = p_node_4;
            element_nodes_9[1] = p_node_1;
            element_nodes_9[2] = p_node_6;
            element_nodes_9[3] = p_node_3;
            Tetrahedra3D4 <NodeType> tetrahedra_9( PointerVector<NodeType>{element_nodes_9} );
            
            std::vector<NodeType::Pointer> element_nodes_10 (4);
            element_nodes_10[0] = p_node_9;
            element_nodes_10[1] = p_node_12;
            element_nodes_10[2] = p_node_11;
            element_nodes_10[3] = p_node_8;
            Tetrahedra3D4 <NodeType> tetrahedra_10( PointerVector<NodeType>{element_nodes_10} );
            
            std::vector<NodeType::Pointer> element_nodes_11 (4);
            element_nodes_11[0] = p_node_3;
            element_nodes_11[1] = p_node_2;
            element_nodes_11[2] = p_node_1;
            element_nodes_11[3] = p_node_6;
            Tetrahedra3D4 <NodeType> tetrahedra_11( PointerVector<NodeType>{element_nodes_11} );
            
            Element::Pointer p_elem_0 = this_model_part.CreateNewElement("Element3D4N", 1, tetrahedra_0, p_elem_prop);
            Element::Pointer p_elem_1 = this_model_part.CreateNewElement("Element3D4N", 2, tetrahedra_1, p_elem_prop);
            Element::Pointer p_elem_2 = this_model_part.CreateNewElement("Element3D4N", 3, tetrahedra_2, p_elem_prop);
            Element::Pointer p_elem_3 = this_model_part.CreateNewElement("Element3D4N", 4, tetrahedra_3, p_elem_prop);
            Element::Pointer p_elem_4 = this_model_part.CreateNewElement("Element3D4N", 5, tetrahedra_4, p_elem_prop);
            Element::Pointer p_elem_5 = this_model_part.CreateNewElement("Element3D4N", 6, tetrahedra_5, p_elem_prop);
            Element::Pointer p_elem_6 = this_model_part.CreateNewElement("Element3D4N", 7, tetrahedra_6, p_elem_prop);
            Element::Pointer p_elem_7 = this_model_part.CreateNewElement("Element3D4N", 8, tetrahedra_7, p_elem_prop);
            Element::Pointer p_elem_8 = this_model_part.CreateNewElement("Element3D4N", 9, tetrahedra_8, p_elem_prop);
            Element::Pointer p_elem_9 = this_model_part.CreateNewElement("Element3D4N", 10, tetrahedra_9, p_elem_prop);
            Element::Pointer p_elem_10 = this_model_part.CreateNewElement("Element3D4N", 11, tetrahedra_10, p_elem_prop);
            Element::Pointer p_elem_11 = this_model_part.CreateNewElement("Element3D4N", 12, tetrahedra_11, p_elem_prop);
            
            // We create the locator
            auto point_locator = BinBasedFastPointLocator<3>(this_model_part);
            point_locator.UpdateSearchDatabase();

            array_1d<double, 3> coordinates(3, 0.0);
            coordinates[0] = 0.5;
            coordinates[1] = 0.5;
            coordinates[2] = 0.5;
            Vector shape_functions;
            Element::Pointer p_element;
            bool is_found = point_locator.FindPointOnMeshSimplified(coordinates, shape_functions, p_element, 1000, 5.0e-2);

            Vector ref_shape_functions(4);
            ref_shape_functions[0] = 0.0;
            ref_shape_functions[1] = 0.5;
            ref_shape_functions[2] = 0.0;
            ref_shape_functions[3] = 0.5;

            const double tolerance = 1.0e-16;
            KRATOS_CHECK(is_found);
            KRATOS_CHECK_EQUAL(p_element->Id(), p_elem_3->Id());
            KRATOS_CHECK_LESS_EQUAL(norm_2(shape_functions - ref_shape_functions), tolerance);

            coordinates[0] = -0.5;
            coordinates[1] = -0.5;
            coordinates[2] = -0.5;
            is_found = point_locator.FindPointOnMeshSimplified(coordinates, shape_functions, p_element, 1000, 5.0e-2);
            KRATOS_CHECK_IS_FALSE(is_found);
        }
        
        /** 
        * Checks the correct work of the binbased fast point locator
        * Test triangle for conditions
        */

        KRATOS_TEST_CASE_IN_SUITE(TestBinBasedFastPointLocator3, KratosCoreFastSuite)
        {
            Model current_model;
            ModelPart& this_model_part = current_model.CreateModelPart("test_model_part",2);
            
            Properties::Pointer p_cond_prop = this_model_part.pGetProperties(0);
            
            auto& process_info = this_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;
            
            // First we create the nodes 
            NodeType::Pointer p_node_1 = this_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_2 = this_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_3 = this_model_part.CreateNewNode(3, 1.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_4 = this_model_part.CreateNewNode(4, 0.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_5 = this_model_part.CreateNewNode(5, 2.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_6 = this_model_part.CreateNewNode(6, 2.0 , 1.0 , 0.0);
            
            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (3);
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            condition_nodes_0[2] = p_node_3;
            Triangle3D3 <NodeType> triangle_0( PointerVector<NodeType>{condition_nodes_0} );
            
            std::vector<NodeType::Pointer> condition_nodes_1 (3);
            condition_nodes_1[0] = p_node_1;
            condition_nodes_1[1] = p_node_3;
            condition_nodes_1[2] = p_node_4;
            Triangle3D3 <NodeType> triangle_1( PointerVector<NodeType>{condition_nodes_1} );
            
            std::vector<NodeType::Pointer> condition_nodes_2 (3);
            condition_nodes_2[0] = p_node_2;
            condition_nodes_2[1] = p_node_5;
            condition_nodes_2[2] = p_node_3;
            Triangle3D3 <NodeType> triangle_2( PointerVector<NodeType>{condition_nodes_2} );
            
            std::vector<NodeType::Pointer> condition_nodes_3 (3);
            condition_nodes_3[0] = p_node_5;
            condition_nodes_3[1] = p_node_6;
            condition_nodes_3[2] = p_node_3;
            Triangle3D3 <NodeType> triangle_3( PointerVector<NodeType>{condition_nodes_3} );
            
            Condition::Pointer p_cond_0 = this_model_part.CreateNewCondition("SurfaceCondition3D3N", 1, triangle_0, p_cond_prop);
            Condition::Pointer p_cond_1 = this_model_part.CreateNewCondition("SurfaceCondition3D3N", 2, triangle_1, p_cond_prop);
            Condition::Pointer p_cond_2 = this_model_part.CreateNewCondition("SurfaceCondition3D3N", 3, triangle_2, p_cond_prop);
            Condition::Pointer p_cond_3 = this_model_part.CreateNewCondition("SurfaceCondition3D3N", 4, triangle_3, p_cond_prop);
            
            // We create the locator
            auto point_locator = BinBasedFastPointLocatorConditions<3>(this_model_part);
            point_locator.UpdateSearchDatabase();

            array_1d<double, 3> coordinates(3, 0.0);
            coordinates[0] = 0.5;
            coordinates[1] = 0.5;
            Vector shape_functions;
            Condition::Pointer p_condition;
            bool is_found = point_locator.FindPointOnMeshSimplified(coordinates, shape_functions, p_condition, 1000, 5.0e-2);

            Vector ref_shape_functions(3);
            ref_shape_functions[0] = 0.5;
            ref_shape_functions[1] = 0.0;
            ref_shape_functions[2] = 0.5;

            const double tolerance = 1.0e-16;
            KRATOS_CHECK(is_found);
            KRATOS_CHECK_EQUAL(p_condition->Id(), p_cond_0->Id());
            KRATOS_CHECK_LESS_EQUAL(norm_2(shape_functions - ref_shape_functions), tolerance);

            coordinates[0] = -0.5;
            coordinates[1] = -0.5;
            is_found = point_locator.FindPointOnMeshSimplified(coordinates, shape_functions, p_condition, 1000, 5.0e-2);
            KRATOS_CHECK_IS_FALSE(is_found);
        }
    } // namespace Testing
}  // namespace Kratos.
