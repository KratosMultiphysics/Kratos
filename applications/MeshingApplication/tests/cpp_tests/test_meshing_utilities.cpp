// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "containers/model.h"
#include "custom_utilities/meshing_utilities.h"

namespace Kratos {
    namespace Testing {

        typedef Node<3> NodeType;
        typedef Geometry<NodeType> GeometryType;

        void CreateDummy2DModelPart(ModelPart& ThisModelPart)
        {
            Properties::Pointer p_elem_prop = ThisModelPart.CreateNewProperties(0);

            // First we create the nodes
            NodeType::Pointer p_node_1 = ThisModelPart.CreateNewNode(1, 0.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_2 = ThisModelPart.CreateNewNode(2, 1.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_3 = ThisModelPart.CreateNewNode(3, 1.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_4 = ThisModelPart.CreateNewNode(4, 0.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_5 = ThisModelPart.CreateNewNode(5, 2.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_6 = ThisModelPart.CreateNewNode(6, 2.0 , 1.0 , 0.0);

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> element_nodes_0 (3);
            element_nodes_0[0] = p_node_1;
            element_nodes_0[1] = p_node_2;
            element_nodes_0[2] = p_node_3;

            std::vector<NodeType::Pointer> element_nodes_1 (3);
            element_nodes_1[0] = p_node_1;
            element_nodes_1[1] = p_node_3;
            element_nodes_1[2] = p_node_4;

            std::vector<NodeType::Pointer> element_nodes_2 (3);
            element_nodes_2[0] = p_node_2;
            element_nodes_2[1] = p_node_5;
            element_nodes_2[2] = p_node_3;

            std::vector<NodeType::Pointer> element_nodes_3 (3);
            element_nodes_3[0] = p_node_5;
            element_nodes_3[1] = p_node_6;
            element_nodes_3[2] = p_node_3;

            Element::Pointer p_elem_0 = ThisModelPart.CreateNewElement("Element2D3N", 1, PointerVector<NodeType>{element_nodes_0}, p_elem_prop);
            Element::Pointer p_elem_1 = ThisModelPart.CreateNewElement("Element2D3N", 2, PointerVector<NodeType>{element_nodes_1}, p_elem_prop);
            Element::Pointer p_elem_2 = ThisModelPart.CreateNewElement("Element2D3N", 3, PointerVector<NodeType>{element_nodes_2}, p_elem_prop);
            Element::Pointer p_elem_3 = ThisModelPart.CreateNewElement("Element2D3N", 4, PointerVector<NodeType>{element_nodes_3}, p_elem_prop);
        }

        void CreateDummy3DModelPart(ModelPart& ThisModelPart)
        {
            Properties::Pointer p_elem_prop = ThisModelPart.CreateNewProperties(0);

            // First we create the nodes
            NodeType::Pointer p_node_1 = ThisModelPart.CreateNewNode(1 , 0.0 , 1.0 , 1.0);
            NodeType::Pointer p_node_2 = ThisModelPart.CreateNewNode(2 , 0.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_3 = ThisModelPart.CreateNewNode(3 , 0.0 , 0.0 , 1.0);
            NodeType::Pointer p_node_4 = ThisModelPart.CreateNewNode(4 , 1.0 , 1.0 , 1.0);
            NodeType::Pointer p_node_5 = ThisModelPart.CreateNewNode(5 , 0.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_6 = ThisModelPart.CreateNewNode(6 , 1.0 , 1.0 , 0.0);

            NodeType::Pointer p_node_7 = ThisModelPart.CreateNewNode(7 , 1.0 , 0.0 , 1.0);
            NodeType::Pointer p_node_8 = ThisModelPart.CreateNewNode(8 , 1.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_9 = ThisModelPart.CreateNewNode(9 , 2.0 , 1.0 , 1.0);
            NodeType::Pointer p_node_10 = ThisModelPart.CreateNewNode(10 , 2.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_11 = ThisModelPart.CreateNewNode(11 , 2.0 , 0.0 , 1.0);
            NodeType::Pointer p_node_12 = ThisModelPart.CreateNewNode(12 , 2.0 , 0.0 , 0.0);

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> element_nodes_0 (4);
            element_nodes_0[0] = p_node_12;
            element_nodes_0[1] = p_node_10;
            element_nodes_0[2] = p_node_8;
            element_nodes_0[3] = p_node_9;

            std::vector<NodeType::Pointer> element_nodes_1 (4);
            element_nodes_1[0] = p_node_4;
            element_nodes_1[1] = p_node_6;
            element_nodes_1[2] = p_node_9;
            element_nodes_1[3] = p_node_7;

            std::vector<NodeType::Pointer> element_nodes_2 (4);
            element_nodes_2[0] = p_node_11;
            element_nodes_2[1] = p_node_7;
            element_nodes_2[2] = p_node_9;
            element_nodes_2[3] = p_node_8;

            std::vector<NodeType::Pointer> element_nodes_3 (4);
            element_nodes_3[0] = p_node_5;
            element_nodes_3[1] = p_node_3;
            element_nodes_3[2] = p_node_8;
            element_nodes_3[3] = p_node_6;

            std::vector<NodeType::Pointer> element_nodes_4 (4);
            element_nodes_4[0] = p_node_4;
            element_nodes_4[1] = p_node_6;
            element_nodes_4[2] = p_node_7;
            element_nodes_4[3] = p_node_3;

            std::vector<NodeType::Pointer> element_nodes_5 (4);
            element_nodes_5[0] = p_node_2;
            element_nodes_5[1] = p_node_3;
            element_nodes_5[2] = p_node_5;
            element_nodes_5[3] = p_node_6;

            std::vector<NodeType::Pointer> element_nodes_6 (4);
            element_nodes_6[0] = p_node_10;
            element_nodes_6[1] = p_node_9;
            element_nodes_6[2] = p_node_6;
            element_nodes_6[3] = p_node_8;

            std::vector<NodeType::Pointer> element_nodes_7 (4);
            element_nodes_7[0] = p_node_7;
            element_nodes_7[1] = p_node_8;
            element_nodes_7[2] = p_node_3;
            element_nodes_7[3] = p_node_6;

            std::vector<NodeType::Pointer> element_nodes_8 (4);
            element_nodes_8[0] = p_node_7;
            element_nodes_8[1] = p_node_8;
            element_nodes_8[2] = p_node_6;
            element_nodes_8[3] = p_node_9;

            std::vector<NodeType::Pointer> element_nodes_9 (4);
            element_nodes_9[0] = p_node_4;
            element_nodes_9[1] = p_node_1;
            element_nodes_9[2] = p_node_6;
            element_nodes_9[3] = p_node_3;

            std::vector<NodeType::Pointer> element_nodes_10 (4);
            element_nodes_10[0] = p_node_9;
            element_nodes_10[1] = p_node_12;
            element_nodes_10[2] = p_node_11;
            element_nodes_10[3] = p_node_8;

            std::vector<NodeType::Pointer> element_nodes_11 (4);
            element_nodes_11[0] = p_node_3;
            element_nodes_11[1] = p_node_2;
            element_nodes_11[2] = p_node_1;
            element_nodes_11[3] = p_node_6;

            Element::Pointer p_elem_0 = ThisModelPart.CreateNewElement("Element3D4N", 1, PointerVector<NodeType>{element_nodes_0}, p_elem_prop);
            Element::Pointer p_elem_1 = ThisModelPart.CreateNewElement("Element3D4N", 2, PointerVector<NodeType>{element_nodes_1}, p_elem_prop);
            Element::Pointer p_elem_2 = ThisModelPart.CreateNewElement("Element3D4N", 3, PointerVector<NodeType>{element_nodes_2}, p_elem_prop);
            Element::Pointer p_elem_3 = ThisModelPart.CreateNewElement("Element3D4N", 4, PointerVector<NodeType>{element_nodes_3}, p_elem_prop);
            Element::Pointer p_elem_4 = ThisModelPart.CreateNewElement("Element3D4N", 5, PointerVector<NodeType>{element_nodes_4}, p_elem_prop);
            Element::Pointer p_elem_5 = ThisModelPart.CreateNewElement("Element3D4N", 6, PointerVector<NodeType>{element_nodes_5}, p_elem_prop);
            Element::Pointer p_elem_6 = ThisModelPart.CreateNewElement("Element3D4N", 7, PointerVector<NodeType>{element_nodes_6}, p_elem_prop);
            Element::Pointer p_elem_7 = ThisModelPart.CreateNewElement("Element3D4N", 8, PointerVector<NodeType>{element_nodes_7}, p_elem_prop);
            Element::Pointer p_elem_8 = ThisModelPart.CreateNewElement("Element3D4N", 9, PointerVector<NodeType>{element_nodes_8}, p_elem_prop);
            Element::Pointer p_elem_9 = ThisModelPart.CreateNewElement("Element3D4N", 10, PointerVector<NodeType>{element_nodes_9}, p_elem_prop);
            Element::Pointer p_elem_10 = ThisModelPart.CreateNewElement("Element3D4N", 11, PointerVector<NodeType>{element_nodes_10}, p_elem_prop);
            Element::Pointer p_elem_11 = ThisModelPart.CreateNewElement("Element3D4N", 12, PointerVector<NodeType>{element_nodes_11}, p_elem_prop);
        }

        /**
        * Checks the correct work of the BlockThresholdSizeElements
        * Test triangle
        */
        KRATOS_TEST_CASE_IN_SUITE(BlockThresholdSizeElements2D, KratosMeshingApplicationFastSuite)
        {
            Model this_model;
            ModelPart& this_model_part = this_model.CreateModelPart("Main", 2);
            ProcessInfo& current_process_info = this_model_part.GetProcessInfo();
            current_process_info[DOMAIN_SIZE] = 2;

            CreateDummy2DModelPart(this_model_part);

            Parameters parameters = Parameters(R"(
            {
                "minimal_size" : 2.0,
                "maximal_size" : 10.0
            })" );

            MeshingUtilities::BlockThresholdSizeElements(this_model_part, parameters);

            for (auto& r_element: this_model_part.Elements()) {
                KRATOS_CHECK(r_element.Is(BLOCKED));
//                 KRATOS_WATCH(r_element.GetValue(ELEMENT_H))
            }
        }

        /**
        * Checks the correct work of the BlockThresholdSizeElements
        * Test tetrahedra
        */
        KRATOS_TEST_CASE_IN_SUITE(BlockThresholdSizeElements3D, KratosMeshingApplicationFastSuite)
        {
            Model this_model;
            ModelPart& this_model_part = this_model.CreateModelPart("Main", 2);
            ProcessInfo& current_process_info = this_model_part.GetProcessInfo();
            current_process_info[DOMAIN_SIZE] = 3;

            CreateDummy3DModelPart(this_model_part);

            Parameters parameters = Parameters(R"(
            {
                "minimal_size" : 2.0,
                "maximal_size" : 10.0
            })" );

            MeshingUtilities::BlockThresholdSizeElements(this_model_part, parameters);

            for (auto& r_element: this_model_part.Elements()) {
                KRATOS_CHECK(r_element.Is(BLOCKED));
//                 KRATOS_WATCH(r_element.GetValue(ELEMENT_H))
            }
        }

    } // namespace Testing
} // namespace Kratos
