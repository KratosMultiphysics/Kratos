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
#include <unordered_set>

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "includes/kratos_flags.h"
#include "includes/kratos_filesystem.h"
#include "utilities/cpp_tests_utilities.h"

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
            ModelPart& r_first_model_part = current_model.CreateModelPart("Main");
            ModelPart& r_first_sub_modelpart_1 = r_first_model_part.CreateSubModelPart("BSubModelPart1");
            ModelPart& r_first_sub_modelpart_2 = r_first_model_part.CreateSubModelPart("ASubModelPart2");
            ModelPart& r_first_sub_modelpart_3 = r_first_model_part.CreateSubModelPart("ZSubModelPart3");
            ModelPart& r_first_sub_modelpart_4 = r_first_model_part.CreateSubModelPart("YSubModelPart4");

            CppTestsUtilities::Create2DGeometry(r_first_model_part, "Element2D3N");

            // Adding nodes to random submodelparts
            r_first_sub_modelpart_1.AddNode(r_first_model_part.pGetNode(1));
            r_first_sub_modelpart_2.AddNode(r_first_model_part.pGetNode(4));
            r_first_sub_modelpart_3.AddNode(r_first_model_part.pGetNode(5));
            r_first_sub_modelpart_4.AddNode(r_first_model_part.pGetNode(6));

            // Adding nodes to random submodelparts
            r_first_sub_modelpart_1.AddElement(r_first_model_part.pGetElement(1));
            r_first_sub_modelpart_2.AddElement(r_first_model_part.pGetElement(4));

//             // Debug
//             KRATOS_WATCH(first_model_part)

            AssignUniqueModelPartCollectionTagUtility collections_utility(r_first_model_part);

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
            second_model_part.AddNodes(r_first_model_part.Nodes().begin(), r_first_model_part.Nodes().end());
            second_model_part.AddElements(r_first_model_part.Elements().begin(), r_first_model_part.Elements().end());

            for (auto & r_nodes_tag : nodes_tags) {
                const IndexSize id = r_nodes_tag.first;
                NodeType::Pointer p_node = second_model_part.pGetNode(id);
                const IndexSize key = r_nodes_tag.second;
                if (key != 0) {// NOTE: key == 0 is the MainModelPart
                    if (collections.find(key) != collections.end()) {
                        for (auto sub_model_part_name : collections[key]) {
                            ModelPart& r_sub_model_part = second_model_part.GetSubModelPart(sub_model_part_name);
                            r_sub_model_part.AddNode(p_node);
                        }
                    }
                }
            }
            for (auto & r_elems_tag : elems_tags) {
                const IndexSize id = r_elems_tag.first;
                ElementType::Pointer p_elem = second_model_part.pGetElement(id);
                const IndexSize key = r_elems_tag.second;
                if (key != 0) {// NOTE: key == 0 is the MainModelPart
                    if (collections.find(key) != collections.end()) {
                        for (auto sub_model_part_name : collections[key]) {
                            ModelPart& r_sub_model_part = second_model_part.GetSubModelPart(sub_model_part_name);
                            r_sub_model_part.AddElement(p_elem);
                        }
                    }
                }
            }

            std::vector<std::string> sub_model_parts_names = r_first_model_part.GetSubModelPartNames();
            for (auto& r_sub_model_part_name : sub_model_parts_names) {
                ModelPart& r_first_sub_model_part = r_first_model_part.GetSubModelPart(r_sub_model_part_name);
                ModelPart& r_second_sub_model_part = second_model_part.GetSubModelPart(r_sub_model_part_name);
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
            ModelPart& r_first_model_part = current_model.CreateModelPart("Main");
            ModelPart& r_first_sub_modelpart_1 = r_first_model_part.CreateSubModelPart("BSubModelPart1");
            ModelPart& r_first_sub_modelpart_1a = r_first_sub_modelpart_1.CreateSubModelPart("SubModelPart1a");
            ModelPart& r_first_sub_modelpart_1b = r_first_sub_modelpart_1.CreateSubModelPart("SubModelPart1b");
            ModelPart& r_first_sub_modelpart_2 = r_first_model_part.CreateSubModelPart("ASubModelPart2");
            ModelPart& r_first_sub_modelpart_3 = r_first_model_part.CreateSubModelPart("ZSubModelPart3");
            ModelPart& r_first_sub_modelpart_4 = r_first_model_part.CreateSubModelPart("YSubModelPart4");

            CppTestsUtilities::Create2DGeometry(r_first_model_part, "Element2D3N");

            // Adding nodes to random submodelparts
            r_first_sub_modelpart_1.AddNode(r_first_model_part.pGetNode(1));
            r_first_sub_modelpart_1a.AddNode(r_first_model_part.pGetNode(2));
            r_first_sub_modelpart_1b.AddNode(r_first_model_part.pGetNode(3));
            r_first_sub_modelpart_2.AddNode(r_first_model_part.pGetNode(4));
            r_first_sub_modelpart_3.AddNode(r_first_model_part.pGetNode(5));
            r_first_sub_modelpart_4.AddNode(r_first_model_part.pGetNode(6));

            // Adding nodes to random submodelparts
            r_first_sub_modelpart_1.AddElement(r_first_model_part.pGetElement(1));
            r_first_sub_modelpart_1a.AddElement(r_first_model_part.pGetElement(2));
            r_first_sub_modelpart_1b.AddElement(r_first_model_part.pGetElement(3));
            r_first_sub_modelpart_2.AddElement(r_first_model_part.pGetElement(4));

            AssignUniqueModelPartCollectionTagUtility collections_utility(r_first_model_part);

            IndexIndexMapType nodes_tags, conds_tags, elems_tags;
            IndexStringMapType collections;
            collections_utility.ComputeTags(nodes_tags, conds_tags, elems_tags, collections);

            // Creating the second model part
            ModelPart& r_second_model_part = current_model.CreateModelPart("SecondMain");
            ModelPart* p_second_sub_modelpart_1 = &r_second_model_part.CreateSubModelPart("BSubModelPart1");
            p_second_sub_modelpart_1->CreateSubModelPart("SubModelPart1a");
            p_second_sub_modelpart_1->CreateSubModelPart("SubModelPart1b");
            r_second_model_part.CreateSubModelPart("ASubModelPart2");
            r_second_model_part.CreateSubModelPart("ZSubModelPart3");
            r_second_model_part.CreateSubModelPart("YSubModelPart4");

            // We add the nodes and elements from the first model part
            r_second_model_part.AddNodes(r_first_model_part.Nodes().begin(), r_first_model_part.Nodes().end());
            r_second_model_part.AddElements(r_first_model_part.Elements().begin(), r_first_model_part.Elements().end());

            for (auto & r_nodes_tag : nodes_tags) {
                const IndexSize id = r_nodes_tag.first;
                NodeType::Pointer p_node = r_second_model_part.pGetNode(id);
                const IndexSize key = r_nodes_tag.second;
                if (key != 0) {// NOTE: key == 0 is the MainModelPart
                    if (collections.find(key) != collections.end()) {
                        for (auto sub_model_part_name : collections[key]) {
                            ModelPart& r_sub_model_part = AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPart(r_second_model_part, sub_model_part_name);
                            r_sub_model_part.AddNode(p_node);
                        }
                    }
                }
            }
            for (auto & r_elems_tag : elems_tags) {
                const IndexSize id = r_elems_tag.first;
                ElementType::Pointer p_elem = r_second_model_part.pGetElement(id);
                const IndexSize key = r_elems_tag.second;
                if (key != 0) {// NOTE: key == 0 is the MainModelPart
                    if (collections.find(key) != collections.end()) {
                        for (auto sub_model_part_name : collections[key]) {
                            ModelPart& r_sub_model_part = AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPart(r_second_model_part, sub_model_part_name);
                            r_sub_model_part.AddElement(p_elem);
                        }
                    }
                }
            }

            std::vector<std::string> sub_model_parts_names = AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPartNames(r_first_model_part);
            for (auto& r_sub_model_part_name : sub_model_parts_names) {
                if (r_sub_model_part_name != r_first_model_part.Name())
                {
                    ModelPart& r_first_sub_model_part = AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPart(r_first_model_part, r_sub_model_part_name);
                    ModelPart& r_second_sub_model_part = AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPart(r_second_model_part, r_sub_model_part_name);
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
            ModelPart& r_model_part = current_model.CreateModelPart("Main");
            ModelPart& r_sub_modelpart_1 = r_model_part.CreateSubModelPart("BSubModelPart1");
            r_sub_modelpart_1.CreateSubModelPart("SubModelPart1a");
            r_sub_modelpart_1.CreateSubModelPart("SubModelPart1b");
            r_model_part.CreateSubModelPart("ASubModelPart2");
            r_model_part.CreateSubModelPart("ZSubModelPart3");
            r_model_part.CreateSubModelPart("YSubModelPart4");

            // auto names = AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPartNames(model_part);
            // for(auto name : names)
            //     KRATOS_WATCH(name);

            KRATOS_CHECK_EQUAL("Main",
                AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPart(r_model_part, "Main").Name());
            KRATOS_CHECK_EQUAL("BSubModelPart1",
                AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPart(r_model_part, "BSubModelPart1").Name());
            KRATOS_CHECK_EQUAL("BSubModelPart1",
                AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPart(r_model_part, "Main.BSubModelPart1").Name());
            KRATOS_CHECK_EQUAL("SubModelPart1b",
                AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPart(r_model_part, "BSubModelPart1.SubModelPart1b").Name());
        }

        /**
        * Checks the correct work of the model parts collections utility static methods
        */
        KRATOS_TEST_CASE_IN_SUITE(ReadTagsFromJsonWriteTagsToJson, KratosCoreFastSuite)
        {
            // Creating the reference model part and the relative submodelparts
            Model current_model;
            ModelPart& r_model_part = current_model.CreateModelPart("Main");
            ModelPart& r_sub_modelpart_1 = r_model_part.CreateSubModelPart("BSubModelPart1");
            r_sub_modelpart_1.CreateSubModelPart("SubModelPart1a");
            r_sub_modelpart_1.CreateSubModelPart("SubModelPart1b");
            r_model_part.CreateSubModelPart("ASubModelPart2");
            r_model_part.CreateSubModelPart("ZSubModelPart3");
            r_model_part.CreateSubModelPart("YSubModelPart4");

            AssignUniqueModelPartCollectionTagUtility collections_utility(r_model_part);

            IndexIndexMapType nodes_tags, conds_tags, elems_tags;
            IndexStringMapType collections;
            collections_utility.ComputeTags(nodes_tags, conds_tags, elems_tags, collections);

            const std::string filename = "test";

            const Parameters param_write = AssignUniqueModelPartCollectionTagUtility::WriteTagsToJson(filename, collections);
            const Parameters param_read = AssignUniqueModelPartCollectionTagUtility::ReadTagsFromJson(filename, collections);

            for (auto itr = param_write.begin(); itr != param_write.end(); ++itr) {
                const std::string& r_name = itr.name();
                KRATOS_CHECK(param_read.Has(r_name));
                const auto& r_write_string_array = itr->GetStringArray();
                const auto& r_read_string_array = param_read[r_name].GetStringArray();
                std::unordered_set<std::string> aux_set;
                for (auto& r_name : r_read_string_array) {
                    aux_set.insert(r_name);
                }
                for (auto& r_name : r_write_string_array) {
                    KRATOS_CHECK(aux_set.find(r_name) != aux_set.end());
                }
            }

            Kratos::filesystem::remove(FilesystemExtensions::JoinPaths({FilesystemExtensions::CurrentWorkingDirectory(), filename + ".json"}));
        }

    } // namespace Testing
}  // namespace Kratos.
