////    |  /           |
////    ' /   __| _` | __|  _ \   __|
////    . \  |   (   | |   (   |\__ `
////   _|\_\_|  \__,_|\__|\___/ ____/
////                   Multi-Physics
////
////  License:         BSD License
////                     Kratos default license: kratos/license.txt
////
////  Main authors:    Tobias Teschemacher
////                   Philipp Bucher
////
//
//// Project includes
//#include "containers/model.h"
//#include "testing/testing.h"
//#include "includes/model_part.h"
//#include "utilities/auxiliar_model_part_utilities.h"
//
//namespace Kratos {
//  namespace Testing {
//
//    typedef Node<3> NodeType;
//
//    void GenerateGenericModelPart(ModelPart& rModelPart)
//    {
//        Properties::Pointer p_elem_prop = rModelPart.CreateNewProperties(0);
//
//        // First we create the nodes
//        NodeType::Pointer p_node_1 = rModelPart.CreateNewNode(1, 0.0 , 0.0 , 0.0);
//        NodeType::Pointer p_node_2 = rModelPart.CreateNewNode(2, 1.0 , 0.0 , 0.0);
//        NodeType::Pointer p_node_3 = rModelPart.CreateNewNode(3, 1.0 , 1.0 , 0.0);
//        NodeType::Pointer p_node_4 = rModelPart.CreateNewNode(4, 0.0 , 1.0 , 0.0);
//        NodeType::Pointer p_node_5 = rModelPart.CreateNewNode(5, 2.0 , 0.0 , 0.0);
//        NodeType::Pointer p_node_6 = rModelPart.CreateNewNode(6, 2.0 , 1.0 , 0.0);
//
//        // Now we create the "conditions"
//        std::vector<NodeType::Pointer> condition_nodes_0 (2);
//        condition_nodes_0[0] = p_node_1;
//        condition_nodes_0[1] = p_node_2;
//
//        std::vector<NodeType::Pointer> condition_nodes_1 (2);
//        condition_nodes_1[0] = p_node_1;
//        condition_nodes_1[1] = p_node_4;
//
//        std::vector<NodeType::Pointer> condition_nodes_2 (2);
//        condition_nodes_2[0] = p_node_2;
//        condition_nodes_2[1] = p_node_5;
//
//        std::vector<NodeType::Pointer> condition_nodes_3 (2);
//        condition_nodes_3[0] = p_node_5;
//        condition_nodes_3[1] = p_node_6;
//
//        Condition::Pointer p_cond_0 = rModelPart.CreateNewCondition("Condition2D2N", 1, PointerVector<NodeType>{condition_nodes_0}, p_elem_prop);
//        Condition::Pointer p_cond_1 = rModelPart.CreateNewCondition("Condition2D2N", 2, PointerVector<NodeType>{condition_nodes_1}, p_elem_prop);
//        Condition::Pointer p_cond_2 = rModelPart.CreateNewCondition("Condition2D2N", 3, PointerVector<NodeType>{condition_nodes_2}, p_elem_prop);
//        Condition::Pointer p_cond_3 = rModelPart.CreateNewCondition("Condition2D2N", 4, PointerVector<NodeType>{condition_nodes_3}, p_elem_prop);
//
//        // Now we create the "elements"
//        std::vector<NodeType::Pointer> element_nodes_0 (3);
//        element_nodes_0[0] = p_node_1;
//        element_nodes_0[1] = p_node_2;
//        element_nodes_0[2] = p_node_3;
//
//        std::vector<NodeType::Pointer> element_nodes_1 (3);
//        element_nodes_1[0] = p_node_1;
//        element_nodes_1[1] = p_node_3;
//        element_nodes_1[2] = p_node_4;
//
//        std::vector<NodeType::Pointer> element_nodes_2 (3);
//        element_nodes_2[0] = p_node_2;
//        element_nodes_2[1] = p_node_5;
//        element_nodes_2[2] = p_node_3;
//
//        std::vector<NodeType::Pointer> element_nodes_3 (3);
//        element_nodes_3[0] = p_node_5;
//        element_nodes_3[1] = p_node_6;
//        element_nodes_3[2] = p_node_3;
//
//        Element::Pointer p_elem_0 = rModelPart.CreateNewElement("Element2D3N", 1, PointerVector<NodeType>{element_nodes_0}, p_elem_prop);
//        Element::Pointer p_elem_1 = rModelPart.CreateNewElement("Element2D3N", 2, PointerVector<NodeType>{element_nodes_1}, p_elem_prop);
//        Element::Pointer p_elem_2 = rModelPart.CreateNewElement("Element2D3N", 3, PointerVector<NodeType>{element_nodes_2}, p_elem_prop);
//        Element::Pointer p_elem_3 = rModelPart.CreateNewElement("Element2D3N", 4, PointerVector<NodeType>{element_nodes_3}, p_elem_prop);
//    }
//
//    KRATOS_TEST_CASE_IN_SUITE(ModelPartGeometryModel, KratosCoreFastSuite)
//    {
//        Model current_model;
//
//        ModelPart& r_model_part = current_model.CreateModelPart("Main");
//
//        r_model_part.CreateSubModelPart("GeometryModel1");
//        r_model_part.CreateSubModelPart("GeometryModel2");
//
//
//
//        KRATOS_CHECK_EQUAL(r_model_part.NumberOfGeometries(), 4);
//        KRATOS_CHECK_EQUAL(r_model_part.GetSubModelPart("GeometryModel1").NumberOfGeometries(), 1);
//        KRATOS_CHECK_EQUAL(r_model_part.GetSubModelPart("GeometryModel2").NumberOfGeometries(), 1);
//    }
//  }  // namespace Testing.
//}  // namespace Kratos.
