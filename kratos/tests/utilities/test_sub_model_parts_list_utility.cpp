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
#include "geometries/triangle_2d_3.h"
#include "testing/testing.h"
#include "includes/kratos_flags.h"

/* Utilities */
#include "utilities/sub_model_parts_list_utility.h"

namespace Kratos
{
    namespace Testing
    {
        typedef Node<3> NodeType;
        typedef std::size_t IndexSize;
        typedef std::unordered_map<IndexSize,int> IndexIntMapType;
        typedef std::unordered_map<int,std::vector<std::string>> IntStringMapType;
        typedef std::map<std::pair<int,int>,int> PairIntMapType;
        typedef std::unordered_map<int,std::vector<ModelPart*>> IntModelPartPtrMapType;

        /**
        * Checks the correct work of the sub modelparts list utility
        */

        KRATOS_TEST_CASE_IN_SUITE(TestSubmodelPartsListUtility, KratosSubModelPartsListUtilityFastSuite)
        {
            // Creating the reference model part and the relative submodelparts non alphabetically ordered
            ModelPart first_model_part("Main");
            ModelPart::Pointer p_first_sub_modelpart_1 = first_model_part.CreateSubModelPart("BSubModelPart1");
            ModelPart::Pointer p_first_sub_modelpart_2 = first_model_part.CreateSubModelPart("ASubModelPart2");
            ModelPart::Pointer p_first_sub_modelpart_3 = first_model_part.CreateSubModelPart("ZSubModelPart3");
            ModelPart::Pointer p_first_sub_modelpart_4 = first_model_part.CreateSubModelPart("YSubModelPart4");

            // Creating the Properties
            Properties::Pointer p_elem_prop = first_model_part.pGetProperties(0);

            // First we create the nodes
            NodeType::Pointer p_node_1 = first_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_2 = first_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_3 = first_model_part.CreateNewNode(3, 1.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_4 = first_model_part.CreateNewNode(4, 0.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_5 = first_model_part.CreateNewNode(5, 2.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_6 = first_model_part.CreateNewNode(6, 2.0 , 1.0 , 0.0);

            // Adding nodes to random submodelparts
            p_first_sub_modelpart_1->AddNode(p_node_1);
            p_first_sub_modelpart_2->AddNode(p_node_4);
            p_first_sub_modelpart_3->AddNode(p_node_5);
            p_first_sub_modelpart_4->AddNode(p_node_6);

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> element_nodes_0 (3);
            element_nodes_0[0] = p_node_1;
            element_nodes_0[1] = p_node_2;
            element_nodes_0[2] = p_node_3;
            Triangle2D3 <NodeType> triangle_0( element_nodes_0 );

            std::vector<NodeType::Pointer> element_nodes_1 (3);
            element_nodes_1[0] = p_node_1;
            element_nodes_1[1] = p_node_3;
            element_nodes_1[2] = p_node_4;
            Triangle2D3 <NodeType> triangle_1( element_nodes_1 );

            std::vector<NodeType::Pointer> element_nodes_2 (3);
            element_nodes_2[0] = p_node_2;
            element_nodes_2[1] = p_node_5;
            element_nodes_2[2] = p_node_3;
            Triangle2D3 <NodeType> triangle_2( element_nodes_2 );

            std::vector<NodeType::Pointer> element_nodes_3 (3);
            element_nodes_3[0] = p_node_5;
            element_nodes_3[1] = p_node_6;
            element_nodes_3[2] = p_node_3;
            Triangle2D3 <NodeType> triangle_3( element_nodes_3 );

            Element::Pointer p_elem_0 = first_model_part.CreateNewElement("Element2D3N", 1, triangle_0, p_elem_prop);
            Element::Pointer p_elem_1 = first_model_part.CreateNewElement("Element2D3N", 2, triangle_1, p_elem_prop);
            Element::Pointer p_elem_2 = first_model_part.CreateNewElement("Element2D3N", 3, triangle_2, p_elem_prop);
            Element::Pointer p_elem_3 = first_model_part.CreateNewElement("Element2D3N", 4, triangle_3, p_elem_prop);

            // Adding nodes to random submodelparts
            p_first_sub_modelpart_1->AddElement(p_elem_0);
            p_first_sub_modelpart_2->AddElement(p_elem_3);

//             // Debug
//             KRATOS_WATCH(first_model_part)

            SubModelPartsListUtility colors_utility(first_model_part);

            IndexIntMapType nodes_colors, cond_colors, elem_colors;
            IntStringMapType colors;
            colors_utility.ComputeSubModelPartsList(nodes_colors, cond_colors, elem_colors, colors);

            // Creating the second model part
            ModelPart second_model_part("Main");
            ModelPart::Pointer p_second_sub_modelpart_1 = second_model_part.CreateSubModelPart("BSubModelPart1");
            ModelPart::Pointer p_second_sub_modelpart_2 = second_model_part.CreateSubModelPart("ASubModelPart2");
            ModelPart::Pointer p_second_sub_modelpart_3 = second_model_part.CreateSubModelPart("ZSubModelPart3");
            ModelPart::Pointer p_second_sub_modelpart_4 = second_model_part.CreateSubModelPart("YSubModelPart4");

            // We add the nodes and elements from the first model part
            second_model_part.AddNodes(first_model_part.Nodes().begin(), first_model_part.Nodes().end());
            second_model_part.AddElements(first_model_part.Elements().begin(), first_model_part.Elements().end());

            for (auto & nodes_color : nodes_colors) {
                const int id = nodes_color.first;
                NodeType::Pointer p_node = second_model_part.pGetNode(id);
                const int key = nodes_color.second;
                if (key != 0) {// NOTE: key == 0 is the MainModelPart
                    if (colors.find(key) != colors.end()) {
                        for (auto sub_model_part_name : colors[key]) {
                            ModelPart& r_sub_model_part = second_model_part.GetSubModelPart(sub_model_part_name);
                            r_sub_model_part.AddNode(p_node);
                        }
                    }
                }
            }
            for (auto & elems_color : elem_colors) {
                const int id = elems_color.first;
                ElementType::Pointer p_elem = second_model_part.pGetElement(id);
                const int key = elems_color.second;
                if (key != 0) {// NOTE: key == 0 is the MainModelPart
                    if (colors.find(key) != colors.end()) {
                        for (auto sub_model_part_name : colors[key]) {
                            ModelPart& r_sub_model_part = second_model_part.GetSubModelPart(sub_model_part_name);
                            r_sub_model_part.AddElement(p_elem);
                        }
                    }
                }
            }

//             // Debug
//             KRATOS_WATCH(second_model_part)

            std::vector<std::string> sub_model_parts_names = first_model_part.GetSubModelPartNames();
            for (auto& sub_model_part_name : sub_model_parts_names) {
                ModelPart& r_first_sub_model_part = first_model_part.GetSubModelPart(sub_model_part_name);
                ModelPart& r_second_sub_model_part = second_model_part.GetSubModelPart(sub_model_part_name);
                KRATOS_CHECK_EQUAL(r_first_sub_model_part.NumberOfNodes(), r_second_sub_model_part.NumberOfNodes());
                KRATOS_CHECK_EQUAL(r_first_sub_model_part.NumberOfElements(), r_second_sub_model_part.NumberOfElements());

                for (std::size_t i = 0; i < r_first_sub_model_part.NumberOfNodes(); i++) {
                    auto it_first_node = r_first_sub_model_part.Nodes().begin() + i;
                    auto it_found_second_node = r_second_sub_model_part.Nodes().find(it_first_node->Id());
                    KRATOS_CHECK_NOT_EQUAL(it_found_second_node, r_second_sub_model_part.NodesEnd());
                }
                for (std::size_t i = 0; i < r_first_sub_model_part.NumberOfElements(); i++) {
                    auto it_first_elem = r_first_sub_model_part.Elements().begin() + i;
                    auto it_found_second_elem = r_second_sub_model_part.Elements().find(it_first_elem->Id());
                    KRATOS_CHECK_NOT_EQUAL(it_found_second_elem, r_second_sub_model_part.ElementsEnd());
                }
            }
        }

        /**
        * Checks the correct work of the modelparts colors utility (with different sublevels of modelparts)
        */

        KRATOS_TEST_CASE_IN_SUITE(TestSubModelPartsListUtilityWithSublevels, KratosSubModelPartsListUtilityFastSuite)
        {
            // Creating the reference model part and the relative submodelparts
            ModelPart first_model_part("Main");
            ModelPart::Pointer p_first_sub_modelpart_1 = first_model_part.CreateSubModelPart("BSubModelPart1");
            ModelPart::Pointer p_first_sub_modelpart_1a = p_first_sub_modelpart_1->CreateSubModelPart("SubModelPart1a");
            ModelPart::Pointer p_first_sub_modelpart_1b = p_first_sub_modelpart_1->CreateSubModelPart("SubModelPart1b");
            ModelPart::Pointer p_first_sub_modelpart_2 = first_model_part.CreateSubModelPart("ASubModelPart2");
            ModelPart::Pointer p_first_sub_modelpart_3 = first_model_part.CreateSubModelPart("ZSubModelPart3");
            ModelPart::Pointer p_first_sub_modelpart_4 = first_model_part.CreateSubModelPart("YSubModelPart4");

            // Creating the Properties
            Properties::Pointer p_elem_prop = first_model_part.pGetProperties(0);

            // First we create the nodes
            NodeType::Pointer p_node_1 = first_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_2 = first_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_3 = first_model_part.CreateNewNode(3, 1.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_4 = first_model_part.CreateNewNode(4, 0.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_5 = first_model_part.CreateNewNode(5, 2.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_6 = first_model_part.CreateNewNode(6, 2.0 , 1.0 , 0.0);

            // Adding nodes to random submodelparts
            p_first_sub_modelpart_1->AddNode(p_node_1);
            p_first_sub_modelpart_1a->AddNode(p_node_2);
            p_first_sub_modelpart_1b->AddNode(p_node_3);
            p_first_sub_modelpart_2->AddNode(p_node_4);
            p_first_sub_modelpart_3->AddNode(p_node_5);
            p_first_sub_modelpart_4->AddNode(p_node_6);

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> element_nodes_0 (3);
            element_nodes_0[0] = p_node_1;
            element_nodes_0[1] = p_node_2;
            element_nodes_0[2] = p_node_3;
            Triangle2D3 <NodeType> triangle_0( element_nodes_0 );

            std::vector<NodeType::Pointer> element_nodes_1 (3);
            element_nodes_1[0] = p_node_1;
            element_nodes_1[1] = p_node_3;
            element_nodes_1[2] = p_node_4;
            Triangle2D3 <NodeType> triangle_1( element_nodes_1 );

            std::vector<NodeType::Pointer> element_nodes_2 (3);
            element_nodes_2[0] = p_node_2;
            element_nodes_2[1] = p_node_5;
            element_nodes_2[2] = p_node_3;
            Triangle2D3 <NodeType> triangle_2( element_nodes_2 );

            std::vector<NodeType::Pointer> element_nodes_3 (3);
            element_nodes_3[0] = p_node_5;
            element_nodes_3[1] = p_node_6;
            element_nodes_3[2] = p_node_3;
            Triangle2D3 <NodeType> triangle_3( element_nodes_3 );

            Element::Pointer p_elem_0 = first_model_part.CreateNewElement("Element2D3N", 1, triangle_0, p_elem_prop);
            Element::Pointer p_elem_1 = first_model_part.CreateNewElement("Element2D3N", 2, triangle_1, p_elem_prop);
            Element::Pointer p_elem_2 = first_model_part.CreateNewElement("Element2D3N", 3, triangle_2, p_elem_prop);
            Element::Pointer p_elem_3 = first_model_part.CreateNewElement("Element2D3N", 4, triangle_3, p_elem_prop);

            // Adding nodes to random submodelparts
            p_first_sub_modelpart_1->AddElement(p_elem_0);
            p_first_sub_modelpart_1a->AddElement(p_elem_1);
            p_first_sub_modelpart_1b->AddElement(p_elem_2);
            p_first_sub_modelpart_2->AddElement(p_elem_3);

//             // Debug
//             KRATOS_WATCH(first_model_part)

            SubModelPartsListUtility colors_utility(first_model_part);

            IndexIntMapType nodes_colors, cond_colors, elem_colors;
            IntStringMapType colors;
            colors_utility.ComputeSubModelPartsList(nodes_colors, cond_colors, elem_colors, colors);

            // Creating the second model part
            ModelPart second_model_part("Main");
            ModelPart::Pointer p_second_sub_modelpart_1 = second_model_part.CreateSubModelPart("BSubModelPart1");
            ModelPart::Pointer p_second_sub_modelpart_1a = p_second_sub_modelpart_1->CreateSubModelPart("SubModelPart1a");
            ModelPart::Pointer p_second_sub_modelpart_1b = p_second_sub_modelpart_1->CreateSubModelPart("SubModelPart1b");
            ModelPart::Pointer p_second_sub_modelpart_2 = second_model_part.CreateSubModelPart("ASubModelPart2");
            ModelPart::Pointer p_second_sub_modelpart_3 = second_model_part.CreateSubModelPart("ZSubModelPart3");
            ModelPart::Pointer p_second_sub_modelpart_4 = second_model_part.CreateSubModelPart("YSubModelPart4");

            // We add the nodes and elements from the first model part
            second_model_part.AddNodes(first_model_part.Nodes().begin(), first_model_part.Nodes().end());
            second_model_part.AddElements(first_model_part.Elements().begin(), first_model_part.Elements().end());

            for (auto & nodes_color : nodes_colors) {
                const int id = nodes_color.first;
                NodeType::Pointer p_node = second_model_part.pGetNode(id);
                const int key = nodes_color.second;
                if (key != 0) {// NOTE: key == 0 is the MainModelPart
                    if (colors.find(key) != colors.end()) {
                        for (auto sub_model_part_name : colors[key]) {
                            ModelPart& r_sub_model_part = SubModelPartsListUtility::GetRecursiveSubModelPart(second_model_part, sub_model_part_name);
                            r_sub_model_part.AddNode(p_node);
                        }
                    }
                }
            }
            for (auto & elems_color : elem_colors) {
                const int id = elems_color.first;
                ElementType::Pointer p_elem = second_model_part.pGetElement(id);
                const int key = elems_color.second;
                if (key != 0) {// NOTE: key == 0 is the MainModelPart
                    if (colors.find(key) != colors.end()) {
                        for (auto sub_model_part_name : colors[key]) {
                            ModelPart& r_sub_model_part = SubModelPartsListUtility::GetRecursiveSubModelPart(second_model_part, sub_model_part_name);
                            r_sub_model_part.AddElement(p_elem);
                        }
                    }
                }
            }

//             // Debug
//             KRATOS_WATCH(second_model_part)

            std::vector<std::string> sub_model_parts_names = SubModelPartsListUtility::GetRecursiveSubModelPartNames(first_model_part);
            for (auto& sub_model_part_name : sub_model_parts_names) {
                ModelPart& r_first_sub_model_part = SubModelPartsListUtility::GetRecursiveSubModelPart(first_model_part, sub_model_part_name);
                ModelPart& r_second_sub_model_part = SubModelPartsListUtility::GetRecursiveSubModelPart(second_model_part, sub_model_part_name);
                KRATOS_CHECK_EQUAL(r_first_sub_model_part.NumberOfNodes(), r_second_sub_model_part.NumberOfNodes());
                KRATOS_CHECK_EQUAL(r_first_sub_model_part.NumberOfElements(), r_second_sub_model_part.NumberOfElements());

                for (std::size_t i = 0; i < r_first_sub_model_part.NumberOfNodes(); i++) {
                    auto it_first_node = r_first_sub_model_part.Nodes().begin() + i;
                    auto it_found_second_node = r_second_sub_model_part.Nodes().find(it_first_node->Id());
                    KRATOS_CHECK_NOT_EQUAL(it_found_second_node, r_second_sub_model_part.NodesEnd());
                }
                for (std::size_t i = 0; i < r_first_sub_model_part.NumberOfElements(); i++) {
                    auto it_first_elem = r_first_sub_model_part.Elements().begin() + i;
                    auto it_found_second_elem = r_second_sub_model_part.Elements().find(it_first_elem->Id());
                    KRATOS_CHECK_NOT_EQUAL(it_found_second_elem, r_second_sub_model_part.ElementsEnd());
                }
            }
        }


        /**
        * Checks the correct work of the modelparts colors
        */

        KRATOS_TEST_CASE_IN_SUITE(TestSubModelPartsListUtilityPointers, KratosSubModelPartsListUtilityFastSuite)
        {
            // Creating the reference model part and the relative submodelparts
            ModelPart model_part("ModelPart");
            ModelPart::Pointer p_sub_modelpart_1 = model_part.CreateSubModelPart("First");
            ModelPart::Pointer p_sub_modelpart_2 = model_part.CreateSubModelPart("Second");
            ModelPart::Pointer p_sub_modelpart_3 = p_sub_modelpart_2->CreateSubModelPart("Third");

            // Creating the nodes
            NodeType::Pointer p_node_1 = model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_2 = model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_3 = model_part.CreateNewNode(3, 2.0 , 0.0 , 0.0);

            // Add the nodes to sub model parts
            p_sub_modelpart_1->AddNode(p_node_1);
            p_sub_modelpart_2->AddNode(p_node_2);
            p_sub_modelpart_3->AddNode(p_node_3);

            // Initialize the utility
            SubModelPartsListUtility colors_utility(model_part);

            // Get the colors
            IndexIntMapType nodes_colors, cond_colors, elem_colors;
            IntStringMapType colors;
            colors_utility.ComputeSubModelPartsList(nodes_colors, cond_colors, elem_colors, colors);

            // Get the pointers
            IntModelPartPtrMapType colors_pointers = SubModelPartsListUtility::GetModelPartColorsPointers(model_part, colors);

            // Check both colormaps have the same size
            KRATOS_CHECK_EQUAL(colors.size(), colors_pointers.size());

            for (auto color : colors)
            {
                // Check both colors have the same size
                std::vector<ModelPart*> p_color = colors_pointers[color.first];
                KRATOS_CHECK_EQUAL(color.second.size(), p_color.size());
                for (IndexSize i = 0; i < color.second.size(); i++)
                {
                    // Check the names
                    KRATOS_CHECK_EQUAL(color.second[i], p_color[i]->Name());
                }
            }
        }


        /**
        * Checks the correct work of the modelparts colors utility (computing the colors intersection)
        */

        KRATOS_TEST_CASE_IN_SUITE(TestSubModelPartsListUtilityIntersections, KratosSubModelPartsListUtilityFastSuite)
        {
            // Creating the reference model part and the relative submodelparts
            ModelPart model_part("Main");
            ModelPart::Pointer p_sub_modelpart_1 = model_part.CreateSubModelPart("BSubModelPart1");
            ModelPart::Pointer p_sub_modelpart_2 = model_part.CreateSubModelPart("ASubModelPart2");
            ModelPart::Pointer p_sub_modelpart_3 = model_part.CreateSubModelPart("ZSubModelPart3");

            // First we create the nodes
            NodeType::Pointer p_node_1 = model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_2 = model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_3 = model_part.CreateNewNode(3, 2.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_4 = model_part.CreateNewNode(4, 0.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_5 = model_part.CreateNewNode(5, 1.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_6 = model_part.CreateNewNode(6, 2.0 , 1.0 , 0.0);

            // Add the nodes to sub model parts
            p_sub_modelpart_1->AddNode(p_node_1);
            p_sub_modelpart_1->AddNode(p_node_4);
            p_sub_modelpart_2->AddNode(p_node_4);
            p_sub_modelpart_2->AddNode(p_node_5);
            p_sub_modelpart_2->AddNode(p_node_6);
            p_sub_modelpart_3->AddNode(p_node_3);
            p_sub_modelpart_3->AddNode(p_node_6);

            // Initialize the utility
            SubModelPartsListUtility colors_utility(model_part);

            // Get the colors
            IndexIntMapType nodes_colors, cond_colors, elem_colors;
            IntStringMapType colors;
            colors_utility.ComputeSubModelPartsList(nodes_colors, cond_colors, elem_colors, colors);

            // Compute the intersections
            PairIntMapType intersections;
            SubModelPartsListUtility::IntersectColors(colors, intersections);
            int key;
            // The intersection gives the main model part
            key = intersections[std::minmax(nodes_colors[p_node_1->Id()], nodes_colors[p_node_3->Id()])];
            KRATOS_CHECK_EQUAL(key, nodes_colors[p_node_2->Id()]);
            KRATOS_CHECK_EQUAL(key, 0);

            // The intersection is a sub model part
            key = intersections[std::minmax(nodes_colors[p_node_4->Id()], nodes_colors[p_node_6->Id()])];
            KRATOS_CHECK_EQUAL(key, nodes_colors[p_node_5->Id()]);

            // The input is included in the intersection
            key = intersections[std::minmax(nodes_colors[p_node_1->Id()], nodes_colors[p_node_4->Id()])];
            KRATOS_CHECK_EQUAL(key, nodes_colors[p_node_1->Id()]);
        }

    } // namespace Testing
}  // namespace Kratos.
