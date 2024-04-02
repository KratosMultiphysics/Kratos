//    |  /           |
//    ' /   _| _` | _|  _ \   _|
//    . \  |   (   | |   (   |\_ `
//   _|\_\_|  \_,_|\_|\__/ __/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//                   Marc Nunez
//

// System includes
#include <filesystem>
#include <unordered_set>

// External includes

// Project includes
#include "containers/model.h"
#include "mpi/testing/mpi_testing.h"
#include "includes/kratos_flags.h"
#include "tests/test_utilities/cpp_tests_utilities.h"
#include "mpi/includes/mpi_communicator.h"

/* Utilities */
#include "utilities/assign_unique_model_part_collection_tag_utility.h"

namespace Kratos::Testing
{
    typedef std::size_t IndexSize;
    typedef std::unordered_map<IndexSize,IndexSize> IndexIndexMapType;
    typedef std::unordered_map<IndexSize,std::vector<std::string>> IndexStringMapType;

    /**
    * Checks the correct work of the sub modelparts list utility
    */
    KRATOS_TEST_CASE_IN_SUITE(AssignMPIUniqueModelPartCollectionTagUtility, KratosMPICoreFastSuite)
    {
        // Creating the reference model part and the relative submodelparts non alphabetically ordered
        Model current_model;
        ModelPart& r_model_part = current_model.CreateModelPart("Main");
        ModelPart& r_sub_modelpart_1 = r_model_part.CreateSubModelPart("BSubModelPart1");
        ModelPart& r_sub_modelpart_2 = r_model_part.CreateSubModelPart("ASubModelPart2");
        ModelPart& r_sub_modelpart_3 = r_model_part.CreateSubModelPart("ZSubModelPart3");
        ModelPart& r_sub_modelpart_4 = r_model_part.CreateSubModelPart("YSubModelPart4");

        Communicator::Pointer pnew_comm = Kratos::make_shared< MPICommunicator >(&r_model_part.GetNodalSolutionStepVariablesList(), Testing::GetDefaultDataCommunicator());
        r_model_part.SetCommunicator(pnew_comm);
        auto& r_data_communicator = r_model_part.GetCommunicator().GetDataCommunicator();
        auto rank = r_data_communicator.Rank();
        auto size = r_data_communicator.Size();

        r_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);
        CppTestsUtilities::Create2DGeometry(r_model_part, "Element2D3N");

        r_sub_modelpart_1.AddNode(r_model_part.pGetNode(1));
        r_sub_modelpart_2.AddNode(r_model_part.pGetNode(4));
        r_sub_modelpart_3.AddNode(r_model_part.pGetNode(5));
        r_sub_modelpart_4.AddNode(r_model_part.pGetNode(6));

        // Adding nodes to random submodelparts
        if (r_data_communicator.IsDistributed() && size > 1) {
            if (rank == 0) {
                r_sub_modelpart_3.AddNode(r_model_part.pGetNode(4));
            }
            else {
                r_sub_modelpart_3.AddNode(r_model_part.pGetNode(6));
                r_sub_modelpart_2.AddNode(r_model_part.pGetNode(5));
                r_sub_modelpart_4.AddNode(r_model_part.pGetNode(5));

            }
        } else {
            r_sub_modelpart_3.AddNode(r_model_part.pGetNode(4));
            r_sub_modelpart_3.AddNode(r_model_part.pGetNode(6));
            r_sub_modelpart_2.AddNode(r_model_part.pGetNode(5));
            r_sub_modelpart_4.AddNode(r_model_part.pGetNode(5));
        }

        // Adding elements to random submodelparts
        r_sub_modelpart_1.AddElement(r_model_part.pGetElement(1));
        r_sub_modelpart_2.AddElement(r_model_part.pGetElement(4));

        AssignUniqueModelPartCollectionTagUtility collections_utility(r_model_part);

        IndexIndexMapType nodes_tags, conds_tags, elems_tags;
        IndexStringMapType collections;
        collections_utility.ComputeTags(nodes_tags, conds_tags, elems_tags, collections);

        const std::string filename = "mpi_test";

        AssignUniqueModelPartCollectionTagUtility::WriteTagsToJson(filename + std::to_string(rank), collections);
        KRATOS_EXPECT_EQ(collections.size(), 8);

        r_data_communicator.Barrier();
        for (int i = 0; i < size; i++) {
            if (i != rank) {
                IndexStringMapType read_collections;
                AssignUniqueModelPartCollectionTagUtility::ReadTagsFromJson(filename + std::to_string(i), read_collections);

                KRATOS_EXPECT_EQ(collections.size(),read_collections.size());

                for (IndexType j = 0; j < read_collections.size(); j++) {
                    KRATOS_EXPECT_EQ(collections[j], read_collections[j]);
                }
            }
        }
        r_data_communicator.Barrier();

        std::filesystem::remove(std::filesystem::current_path() / (filename + std::to_string(rank) + ".json"));
    }
} // namespace Kratos::Testing
