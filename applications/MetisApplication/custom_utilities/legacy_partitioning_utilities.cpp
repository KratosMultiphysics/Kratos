//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//                   Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "legacy_partitioning_utilities.h"


namespace Kratos {

void LegacyPartitioningUtilities::CalculateDomainsGraph(
    IO::GraphType& rDomainsGraph,
    SizeType NumberOfElements,
    IO::ConnectivitiesContainerType& ElementsConnectivities,
    PartitionIndicesType const& NPart,
    PartitionIndicesType const&  EPart )
{
    for(SizeType i_element = 0 ; i_element < NumberOfElements ; i_element++)
        for(std::vector<std::size_t>::iterator i_node = ElementsConnectivities[i_element].begin() ;
                i_node != ElementsConnectivities[i_element].end() ; i_node++)
        {
            SizeType node_rank = NPart[*i_node-1];
            SizeType element_rank = EPart[i_element];
            if(node_rank != element_rank)
            {
                rDomainsGraph(node_rank, element_rank) = 1;
                rDomainsGraph(element_rank, node_rank) = 1;
            }
        }
}

void LegacyPartitioningUtilities::DividingNodes(
    IO::PartitionIndicesContainerType& rNodesAllPartitions,
    IO::ConnectivitiesContainerType& rElementsConnectivities,
    IO::ConnectivitiesContainerType& rConditionsConnectivities,
    const PartitionIndicesType& rNodesPartitions,
    const PartitionIndicesType& rElementsPartitions,
    const PartitionIndicesType& rConditionsPartitions
    )
{
    const SizeType number_of_nodes = rNodesPartitions.size();
    const SizeType number_of_elements = rElementsPartitions.size();
    const SizeType number_of_conditions = rConditionsPartitions.size();

    rNodesAllPartitions.resize(number_of_nodes);

    for(IndexType i_element = 0 ; i_element < number_of_elements ; i_element++) {
        const int element_partition = rElementsPartitions[i_element];

        // For each element in the model loop over its connectivities
        auto& r_element_connectivity = rElementsConnectivities[i_element];
        for(auto it_node = r_element_connectivity.begin(); it_node != r_element_connectivity.end() ; it_node++) {
            // Get global id. We assume that node ids are began with one
            const int my_gid = *it_node-1;

            // Get the partition index for the node i am interested in
            const int node_partition = rNodesPartitions[my_gid];

            // Adding the partition of the element to its nodes
            if(element_partition != node_partition) { // We will add the node_partition once afterward
                rNodesAllPartitions[my_gid].push_back(element_partition);
            }
        }
    }

    for(IndexType i_condition = 0 ; i_condition < number_of_conditions ; i_condition++) {
        const int condition_partition = rConditionsPartitions[i_condition];

        // For each element in the model loop over its connectivities
        auto& r_condition_connectivity = rConditionsConnectivities[i_condition];
        for(auto it_node = r_condition_connectivity.begin(); it_node != r_condition_connectivity.end() ; it_node++) {
            // Get global id. We assume that node ids are began with one
            const int my_gid = *it_node-1;

            // Get the partition index for the node i am interested in
            const int node_partition = rNodesPartitions[my_gid];

            // Adding the partition of the element to its nodes
            if (condition_partition != node_partition)  { // we will add the node_partition once afterward
                rNodesAllPartitions[my_gid].push_back(condition_partition);
            }
        }
    }

    // Adding the nodes partition to their array of partitions and clear the repeated ones
    for(IndexType i_node = 0 ; i_node < number_of_nodes ; i_node++) {
        IO::PartitionIndicesContainerType::value_type& node_partitions = rNodesAllPartitions[i_node];
        node_partitions.push_back(rNodesPartitions[i_node]);

        std::sort(node_partitions.begin(), node_partitions.end()); // TODO: Add parallel sort from STL lib
        IO::PartitionIndicesContainerType::value_type::iterator new_end = std::unique(node_partitions.begin(), node_partitions.end());
        node_partitions.resize(new_end - node_partitions.begin());
    }
}

void LegacyPartitioningUtilities::DividingGeometries(
    IO::PartitionIndicesContainerType& rGeometriesAllPartitions,
    const PartitionIndicesType& rGeometriesPartitions
    )
{
    const SizeType number_of_geometries = rGeometriesPartitions.size();

    rGeometriesAllPartitions.resize(number_of_geometries);

    // Adding the geometry partition to their array of partitions
    for(IndexType i_geometry = 0 ; i_geometry < number_of_geometries ; i_geometry++) {
        rGeometriesAllPartitions[i_geometry].push_back(rGeometriesPartitions[i_geometry]);
    }
}

void LegacyPartitioningUtilities::DividingElements(
    IO::PartitionIndicesContainerType& rElementsAllPartitions,
    PartitionIndicesType const& ElementsPartitions)
{
    SizeType number_of_elements = ElementsPartitions.size();

    rElementsAllPartitions.resize(number_of_elements);

    // adding the elements partition to their array of partitions
    for(SizeType i_element = 0 ; i_element < number_of_elements ; i_element++)
    {
        rElementsAllPartitions[i_element].push_back(ElementsPartitions[i_element]);
    }
}

void LegacyPartitioningUtilities::DividingConditions(
    IO::PartitionIndicesContainerType& rConditionsAllPartitions,
    PartitionIndicesType const& ConditionsPartitions)
{
    SizeType number_of_conditions = ConditionsPartitions.size();

    rConditionsAllPartitions.resize(number_of_conditions);

    // adding the condition partition to their array of partitions
    for(SizeType i_condition = 0 ; i_condition < number_of_conditions ; i_condition++)
    {
        rConditionsAllPartitions[i_condition].push_back(ConditionsPartitions[i_condition]);
    }
}

void LegacyPartitioningUtilities::DividingMasterSlaveConstraints(
    IO::PartitionIndicesContainerType& rMasterSlaveConstraintsAllPartitions,
    const PartitionIndicesType& MasterSlaveConstraintsPartitions
    )
{
    const SizeType number_of_constraints = MasterSlaveConstraintsPartitions.size();

    rMasterSlaveConstraintsAllPartitions.resize(number_of_constraints);

    // Adding the constraint partition to their array of partitions
    for(IndexType i_constraint = 0 ; i_constraint < number_of_constraints ; i_constraint++) {
        rMasterSlaveConstraintsAllPartitions[i_constraint].push_back(MasterSlaveConstraintsPartitions[i_constraint]);
    }
}

void LegacyPartitioningUtilities::ConvertKratosToCSRFormat(
    IO::ConnectivitiesContainerType& KratosFormatNodeConnectivities,
    idxtype** NodeIndices,
    idxtype** NodeConnectivities)
{
    idxtype*& rNodeIndices = *NodeIndices;
    idxtype*& rNodeConnectivities = *NodeConnectivities;

    SizeType num_entries = 0;
    for (IO::ConnectivitiesContainerType::iterator it = KratosFormatNodeConnectivities.begin(); it != KratosFormatNodeConnectivities.end(); ++it)
    {
        num_entries += it->size();
    }

    SizeType num_nodes = KratosFormatNodeConnectivities.size();
    rNodeIndices = new idxtype[num_nodes+1];
    rNodeIndices[0] = 0;
    rNodeConnectivities = new idxtype[num_entries];

    SizeType i = 0;
    SizeType aux_index = 0;

    for (IO::ConnectivitiesContainerType::iterator it = KratosFormatNodeConnectivities.begin(); it != KratosFormatNodeConnectivities.end(); ++it)
    {
        for (std::vector<SizeType>::iterator entry_it = it->begin(); entry_it != it->end(); entry_it++)
            rNodeConnectivities[aux_index++] = (*entry_it - 1); // substract 1 to make Ids start from 0
        rNodeIndices[++i] = aux_index;
    }

}


} // namespace Kratos
