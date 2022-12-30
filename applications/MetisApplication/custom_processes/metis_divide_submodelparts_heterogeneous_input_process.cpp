//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Daniel Diez
//



// System includes

// External includes

// Project includes
#include "custom_processes/metis_divide_submodelparts_heterogeneous_input_process.h"
#include "custom_utilities/legacy_partitioning_utilities.h" // TODO remove

namespace Kratos {


void MetisDivideSubModelPartsHeterogeneousInputProcess::GetNodesPartitions(
    std::vector<idxtype> &rNodePartition,
    SizeType &rNumNodes) {

        KRATOS_TRY

        SizeType num_nodes_in_mesh = BaseType::mrIO.ReadNodesNumber();
        rNodePartition.resize(num_nodes_in_mesh);

        for (auto submodelpart_name : mSettings["sub_model_part_list"].GetStringArray()) {
            // Read nodal graph from input
            IO::ConnectivitiesContainerType kratos_format_node_connectivities;
            std::unordered_set<SizeType> submodelpart_elements_ids;
            std::unordered_set<SizeType> submodelpart_conditions_ids;
            std::vector<idxtype> submodelpart_partition_nodes;

            BaseType::mrIO.ReadSubModelPartElementsAndConditionsIds(
                submodelpart_name,
                submodelpart_elements_ids,
                submodelpart_conditions_ids);

            rNumNodes = BaseType::mrIO.ReadNodalGraphFromEntitiesList(
                kratos_format_node_connectivities,
                submodelpart_elements_ids,
                submodelpart_conditions_ids);

            SizeType id = 0;
            IO::ConnectivitiesContainerType reduced_kratos_format_node_connectivities;
            std::unordered_map<SizeType, SizeType> node_ids_map;

            for (SizeType i = 0; i < kratos_format_node_connectivities.size(); i++) {
                if (kratos_format_node_connectivities[i].size() > 0) {
                    node_ids_map.insert(std::map<SizeType, SizeType>::value_type(i, id++));
                }
            }
            reduced_kratos_format_node_connectivities.resize(node_ids_map.size());
            for (SizeType i = 0; i < kratos_format_node_connectivities.size(); i++) {
                if (kratos_format_node_connectivities[i].size() > 0) {
                    SizeType position = node_ids_map[i];
                    std::vector<SizeType> node_connectivities;
                    for (SizeType node_connectivity : kratos_format_node_connectivities[i]) {
                        SizeType new_connec = node_ids_map[node_connectivity-1];
                        node_connectivities.push_back(new_connec+1);
                    }
                    reduced_kratos_format_node_connectivities[position] = node_connectivities;
                }
            }

            // Write connectivity data in CSR format
            idxtype* node_indices = 0;
            idxtype* node_connectivities = 0;


            LegacyPartitioningUtilities::ConvertKratosToCSRFormat(reduced_kratos_format_node_connectivities, &node_indices, &node_connectivities);

            PartitionNodes(reduced_kratos_format_node_connectivities.size(), node_indices, node_connectivities, submodelpart_partition_nodes);


            for (auto map_element : node_ids_map) {
                rNodePartition[map_element.first] = submodelpart_partition_nodes[map_element.second];
            }

            // Free some memory we no longer need
            delete [] node_indices;
            delete [] node_connectivities;
        }

        mNumNodes = rNodePartition.size();

        KRATOS_CATCH("")

    }


///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.
