//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "mpi/tests/test_utilities/mpi_cpp_test_utilities.h"
#include "mpi/utilities/parallel_fill_communicator.h"
#include "mpi/utilities/gather_modelpart_utility.h"

namespace Kratos::Testing 
{

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(GatherModelPartUtilityGatherEntitiesFromOtherPartitions, KratosMPICoreFastSuite)
{
    // The model part
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");

    // The data communicator
    const DataCommunicator& r_data_communicator = Testing::GetDefaultDataCommunicator();

    // Fill the model part
    MPICppTestUtilities::GenerateDistributedBarStructure(r_model_part, r_data_communicator);

    // The global pointer commmunicator for nodes
    GlobalPointersVector<Node> global_pointers;
    auto p_global_pointer_communicator = MPICppTestUtilities::SynchronizeNodes(r_model_part, global_pointers, r_data_communicator);

    // Define the indices vector
    const std::size_t number_of_gp = global_pointers.size();
    std::vector<std::size_t> indices(number_of_gp);
    std::vector<std::size_t> partition_index(number_of_gp);

    // Call Apply to get the proxy
    auto proxy_id = p_global_pointer_communicator->Apply([](GlobalPointer<Node>& rGP) -> std::size_t {
        return rGP->Id();
    });
    auto proxy_partition = p_global_pointer_communicator->Apply([](GlobalPointer<Node>& rGP) -> std::size_t {
        return rGP->FastGetSolutionStepValue(PARTITION_INDEX);
    });

    // Get the indices
    for(std::size_t i=0; i<number_of_gp; ++i) {
        auto& r_gp = global_pointers(i);
        indices[i] = proxy_id.Get(r_gp);
        partition_index[i] = proxy_partition.Get(r_gp);
    }

    // Indices to bring
    std::unordered_set<std::size_t> set_indices_to_bring_just_ids;
    std::vector<std::pair<std::size_t, std::size_t>> set_indices_to_bring;
    for (std::size_t  i = 0; i < number_of_gp; ++i) {
        if (!r_model_part.HasNode(indices[i])) {
            if (set_indices_to_bring_just_ids.find(indices[i]) == set_indices_to_bring_just_ids.end()) {
                set_indices_to_bring_just_ids.insert(indices[i]);
                set_indices_to_bring.push_back({indices[i], partition_index[i]});
            }
        }
    }
    std::vector<std::size_t> partition_origin;
    partition_origin.reserve(set_indices_to_bring.size());
    std::vector<std::size_t> indices_to_bring;
    indices_to_bring.reserve(set_indices_to_bring.size());
    for (auto& r_pair : set_indices_to_bring) {
        indices_to_bring.push_back(r_pair.first);
        partition_origin.push_back(r_pair.second);
    }

    // Generate map
    std::map<int, std::vector<std::size_t>> nodes_to_bring, elements_to_bring, conditions_to_bring;
    for (std::size_t i = 0; i < indices_to_bring.size(); ++i) {
        auto it_found = nodes_to_bring.find(partition_origin[i]);
        if (it_found != nodes_to_bring.end()) {
            it_found->second.push_back(indices_to_bring[i]);
        } else {
            std::vector<std::size_t> minor_vector(1, indices_to_bring[i]);
            nodes_to_bring.insert({partition_origin[i], minor_vector}); 
        }
    }

    // Bring entities
    GatherModelPartUtility::GatherEntitiesFromOtherPartitions(r_model_part, nodes_to_bring,elements_to_bring, conditions_to_bring);

    // Check the number of nodes (all partitions have 11 nodes)
    KRATOS_EXPECT_EQ(r_model_part.NumberOfNodes(), 11);
}

} // namespace Kratos::Testing