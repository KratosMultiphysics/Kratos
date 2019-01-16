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
#include "testing/testing.h"
#include "includes/kratos_flags.h"

/* Utilities */
#include "utilities/assign_unique_model_part_collection_tag_utility.h"

namespace Kratos
{
    namespace Testing
    {
        typedef Node<3> NodeType;
        typedef std::size_t IndexSize;
        typedef std::unordered_map<IndexSize,IndexSize> IndexIndexMapType;
        typedef std::unordered_map<IndexSize,std::vector<std::string>> IndexStringMapType;

        /**
        * Checks the correct work of the sub modelparts list utility
        */

        KRATOS_TEST_CASE_IN_SUITE(AssignUniqueModelPartCollectionTagUtility, KratosCoreFastSuite)
        {
            
            // Creating the reference model part and the relative submodelparts non alphabetically ordered
            Model current_model;
            ModelPart& first_model_part = current_model.CreateModelPart("Main");
            ModelPart& first_sub_modelpart_1 = first_model_part.CreateSubModelPart("BSubModelPart1");
            ModelPart& first_sub_modelpart_2 = first_model_part.CreateSubModelPart("ASubModelPart2");
            ModelPart& first_sub_modelpart_3 = first_model_part.CreateSubModelPart("ZSubModelPart3");
            ModelPart& first_sub_modelpart_4 = first_model_part.CreateSubModelPart("YSubModelPart4");

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
            first_sub_modelpart_1.AddNode(p_node_1);
            first_sub_modelpart_2.AddNode(p_node_4);
            first_sub_modelpart_3.AddNode(p_node_5);
            first_sub_modelpart_4.AddNode(p_node_6);

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

            Element::Pointer p_elem_0 = first_model_part.CreateNewElement("Element2D3N", 1, triangle_0, p_elem_prop);
            Element::Pointer p_elem_1 = first_model_part.CreateNewElement("Element2D3N", 2, triangle_1, p_elem_prop);
            Element::Pointer p_elem_2 = first_model_part.CreateNewElement("Element2D3N", 3, triangle_2, p_elem_prop);
            Element::Pointer p_elem_3 = first_model_part.CreateNewElement("Element2D3N", 4, triangle_3, p_elem_prop);

            // Adding nodes to random submodelparts
            first_sub_modelpart_1.AddElement(p_elem_0);
            first_sub_modelpart_2.AddElement(p_elem_3);

//             // Debug
//             KRATOS_WATCH(first_model_part)

            AssignUniqueModelPartCollectionTagUtility collections_utility(first_model_part);

            IndexIndexMapType nodes_tags, conds_tags, elems_tags;
            IndexStringMapType collections;
            collections_utility.ComputeTags(nodes_tags, conds_tags, elems_tags, collections);

            // Creating the second model part
            ModelPart& second_model_part = current_model.CreateModelPart("SecondMain");
            second_model_part.CreateSubModelPart("BSubModelPart1");
            second_model_part.CreateSubModelPart("ASubModelPart2");
            second_model_part.CreateSubModelPart("ZSubModelPart3");
            second_model_part.CreateSubModelPart("YSubModelPart4");

            // We add the nodes and elements from the first model part
            second_model_part.AddNodes(first_model_part.Nodes().begin(), first_model_part.Nodes().end());
            second_model_part.AddElements(first_model_part.Elements().begin(), first_model_part.Elements().end());

            for (auto & nodes_tag : nodes_tags) {
                const IndexSize id = nodes_tag.first;
                NodeType::Pointer p_node = second_model_part.pGetNode(id);
                const IndexSize key = nodes_tag.second;
                if (key != 0) {// NOTE: key == 0 is the MainModelPart
                    if (collections.find(key) != collections.end()) {
                        for (auto sub_model_part_name : collections[key]) {
                            ModelPart& r_sub_model_part = second_model_part.GetSubModelPart(sub_model_part_name);
                            r_sub_model_part.AddNode(p_node);
                        }
                    }
                }
            }
            for (auto & elems_tag : elems_tags) {
                const IndexSize id = elems_tag.first;
                ElementType::Pointer p_elem = second_model_part.pGetElement(id);
                const IndexSize key = elems_tag.second;
                if (key != 0) {// NOTE: key == 0 is the MainModelPart
                    if (collections.find(key) != collections.end()) {
                        for (auto sub_model_part_name : collections[key]) {
                            ModelPart& r_sub_model_part = second_model_part.GetSubModelPart(sub_model_part_name);
                            r_sub_model_part.AddElement(p_elem);
                        }
                    }
                }
            }


            std::vector<std::string> sub_model_parts_names = first_model_part.GetSubModelPartNames();
            for (auto& sub_model_part_name : sub_model_parts_names) {
                ModelPart& r_first_sub_model_part = first_model_part.GetSubModelPart(sub_model_part_name);
                ModelPart& r_second_sub_model_part = second_model_part.GetSubModelPart(sub_model_part_name);
                KRATOS_CHECK_EQUAL(r_first_sub_model_part.NumberOfNodes(), r_second_sub_model_part.NumberOfNodes());
                KRATOS_CHECK_EQUAL(r_first_sub_model_part.NumberOfElements(), r_second_sub_model_part.NumberOfElements());

                for (IndexType i = 0; i < r_first_sub_model_part.NumberOfNodes(); i++) {
                    auto it_first_node = r_first_sub_model_part.Nodes().begin() + i;
                    auto it_found_second_node = r_second_sub_model_part.Nodes().find(it_first_node->Id());
                    KRATOS_CHECK_NOT_EQUAL(it_found_second_node, r_second_sub_model_part.NodesEnd());
                }
                for (IndexType i = 0; i < r_first_sub_model_part.NumberOfElements(); i++) {
                    auto it_first_elem = r_first_sub_model_part.Elements().begin() + i;
                    auto it_found_second_elem = r_second_sub_model_part.Elements().find(it_first_elem->Id());
                    KRATOS_CHECK_NOT_EQUAL(it_found_second_elem, r_second_sub_model_part.ElementsEnd());
                }
            }
        }


        /**
        * Checks the correct work of the model parts collections utility (with different sublevels of modelparts)
        */
        KRATOS_TEST_CASE_IN_SUITE(AssignUniqueModelPartCollectionTagUtilityWithSubLevels, KratosCoreFastSuite)
        {
            // Creating the reference model part and the relative submodelparts
            Model current_model;
            ModelPart& first_model_part = current_model.CreateModelPart("Main");
            ModelPart& first_sub_modelpart_1 = first_model_part.CreateSubModelPart("BSubModelPart1");
            ModelPart& first_sub_modelpart_1a = first_sub_modelpart_1.CreateSubModelPart("SubModelPart1a");
            ModelPart& first_sub_modelpart_1b = first_sub_modelpart_1.CreateSubModelPart("SubModelPart1b");
            ModelPart& first_sub_modelpart_2 = first_model_part.CreateSubModelPart("ASubModelPart2");
            ModelPart& first_sub_modelpart_3 = first_model_part.CreateSubModelPart("ZSubModelPart3");
            ModelPart& first_sub_modelpart_4 = first_model_part.CreateSubModelPart("YSubModelPart4");

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
            first_sub_modelpart_1.AddNode(p_node_1);
            first_sub_modelpart_1a.AddNode(p_node_2);
            first_sub_modelpart_1b.AddNode(p_node_3);
            first_sub_modelpart_2.AddNode(p_node_4);
            first_sub_modelpart_3.AddNode(p_node_5);
            first_sub_modelpart_4.AddNode(p_node_6);

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

            Element::Pointer p_elem_0 = first_model_part.CreateNewElement("Element2D3N", 1, triangle_0, p_elem_prop);
            Element::Pointer p_elem_1 = first_model_part.CreateNewElement("Element2D3N", 2, triangle_1, p_elem_prop);
            Element::Pointer p_elem_2 = first_model_part.CreateNewElement("Element2D3N", 3, triangle_2, p_elem_prop);
            Element::Pointer p_elem_3 = first_model_part.CreateNewElement("Element2D3N", 4, triangle_3, p_elem_prop);

            // Adding nodes to random submodelparts
            first_sub_modelpart_1.AddElement(p_elem_0);
            first_sub_modelpart_1a.AddElement(p_elem_1);
            first_sub_modelpart_1b.AddElement(p_elem_2);
            first_sub_modelpart_2.AddElement(p_elem_3);

            AssignUniqueModelPartCollectionTagUtility collections_utility(first_model_part);

            IndexIndexMapType nodes_tags, conds_tags, elems_tags;
            IndexStringMapType collections;
            collections_utility.ComputeTags(nodes_tags, conds_tags, elems_tags, collections);

            // Creating the second model part
            ModelPart& second_model_part = current_model.CreateModelPart("SecondMain");
            ModelPart* p_second_sub_modelpart_1 = &second_model_part.CreateSubModelPart("BSubModelPart1");
            p_second_sub_modelpart_1->CreateSubModelPart("SubModelPart1a");
            p_second_sub_modelpart_1->CreateSubModelPart("SubModelPart1b");
            second_model_part.CreateSubModelPart("ASubModelPart2");
            second_model_part.CreateSubModelPart("ZSubModelPart3");
            second_model_part.CreateSubModelPart("YSubModelPart4");

            // We add the nodes and elements from the first model part
            second_model_part.AddNodes(first_model_part.Nodes().begin(), first_model_part.Nodes().end());
            second_model_part.AddElements(first_model_part.Elements().begin(), first_model_part.Elements().end());

            for (auto & nodes_tag : nodes_tags) {
                const IndexSize id = nodes_tag.first;
                NodeType::Pointer p_node = second_model_part.pGetNode(id);
                const IndexSize key = nodes_tag.second;
                if (key != 0) {// NOTE: key == 0 is the MainModelPart
                    if (collections.find(key) != collections.end()) {
                        for (auto sub_model_part_name : collections[key]) {
                            ModelPart& r_sub_model_part = AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPart(second_model_part, sub_model_part_name);
                            r_sub_model_part.AddNode(p_node);
                        }
                    }
                }
            }
            for (auto & elems_tag : elems_tags) {
                const IndexSize id = elems_tag.first;
                ElementType::Pointer p_elem = second_model_part.pGetElement(id);
                const IndexSize key = elems_tag.second;
                if (key != 0) {// NOTE: key == 0 is the MainModelPart
                    if (collections.find(key) != collections.end()) {
                        for (auto sub_model_part_name : collections[key]) {
                            ModelPart& r_sub_model_part = AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPart(second_model_part, sub_model_part_name);
                            r_sub_model_part.AddElement(p_elem);
                        }
                    }
                }
            }

            std::vector<std::string> sub_model_parts_names = AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPartNames(first_model_part);
            for (auto& sub_model_part_name : sub_model_parts_names) {
                if (sub_model_part_name != first_model_part.Name())
                {
                    ModelPart& r_first_sub_model_part = AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPart(first_model_part, sub_model_part_name);
                    ModelPart& r_second_sub_model_part = AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPart(second_model_part, sub_model_part_name);
                    KRATOS_CHECK_EQUAL(r_first_sub_model_part.NumberOfNodes(), r_second_sub_model_part.NumberOfNodes());
                    KRATOS_CHECK_EQUAL(r_first_sub_model_part.NumberOfElements(), r_second_sub_model_part.NumberOfElements());

                    for (IndexType i = 0; i < r_first_sub_model_part.NumberOfNodes(); i++) {
                        auto it_first_node = r_first_sub_model_part.Nodes().begin() + i;
                        auto it_found_second_node = r_second_sub_model_part.Nodes().find(it_first_node->Id());
                        KRATOS_CHECK_NOT_EQUAL(it_found_second_node, r_second_sub_model_part.NodesEnd());
                    }
                    for (IndexType i = 0; i < r_first_sub_model_part.NumberOfElements(); i++) {
                        auto it_first_elem = r_first_sub_model_part.Elements().begin() + i;
                        auto it_found_second_elem = r_second_sub_model_part.Elements().find(it_first_elem->Id());
                        KRATOS_CHECK_NOT_EQUAL(it_found_second_elem, r_second_sub_model_part.ElementsEnd());
                    }
                }
            }
        }


        /**
        * Checks the correct work of the model parts collections utility static methods
        */
        KRATOS_TEST_CASE_IN_SUITE(GetRecursiveSubModelPart, KratosCoreFastSuite)
        {
            // Creating the reference model part and the relative submodelparts
            Model current_model;
            ModelPart& model_part = current_model.CreateModelPart("Main");
            ModelPart& sub_modelpart_1 = model_part.CreateSubModelPart("BSubModelPart1");
            sub_modelpart_1.CreateSubModelPart("SubModelPart1a");
            sub_modelpart_1.CreateSubModelPart("SubModelPart1b");
            model_part.CreateSubModelPart("ASubModelPart2");
            model_part.CreateSubModelPart("ZSubModelPart3");
            model_part.CreateSubModelPart("YSubModelPart4");

            // auto names = AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPartNames(model_part);
            // for(auto name : names)
            //     KRATOS_WATCH(name);

            KRATOS_CHECK_EQUAL("Main",
                AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPart(model_part, "Main").Name());
            KRATOS_CHECK_EQUAL("BSubModelPart1",
                AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPart(model_part, "BSubModelPart1").Name());
            KRATOS_CHECK_EQUAL("BSubModelPart1",
                AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPart(model_part, "Main.BSubModelPart1").Name());
            KRATOS_CHECK_EQUAL("SubModelPart1b",
                AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPart(model_part, "BSubModelPart1.SubModelPart1b").Name());
        }

    } // namespace Testing
}  // namespace Kratos.
