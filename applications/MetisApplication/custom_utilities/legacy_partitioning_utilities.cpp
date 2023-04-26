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
    IO::ConnectivitiesContainerType& ElementsConnectivities,
    IO::ConnectivitiesContainerType& ConditionsConnectivities,
    PartitionIndicesType const& NodesPartitions,
    PartitionIndicesType const& ElementsPartitions,
    PartitionIndicesType const& ConditionsPartitions)
{
    SizeType number_of_nodes = NodesPartitions.size();
    SizeType number_of_elements = ElementsPartitions.size();
    SizeType number_of_conditions = ConditionsPartitions.size();

    rNodesAllPartitions.resize(number_of_nodes);

    for(SizeType i_element = 0 ; i_element < number_of_elements ; i_element++)
    {
        const int element_partition = ElementsPartitions[i_element];

        //for each element in the model loop over its connectivities
        for(IO::ConnectivitiesContainerType::value_type::iterator i_node = ElementsConnectivities[i_element].begin() ;
                i_node != ElementsConnectivities[i_element].end() ; i_node++)
        {
            //get global id. We assume that node ids are began with one
            const int my_gid = *i_node-1;

            //get the partition index for the node i am interested in
            const int node_partition = NodesPartitions[my_gid];

            // adding the partition of the element to its nodes
            if(element_partition != node_partition) // we will add the node_partition once afterward
                rNodesAllPartitions[my_gid].push_back(element_partition);
        }
    }

    for(SizeType i_condition = 0 ; i_condition < number_of_conditions ; i_condition++)
    {
        const int condition_partition = ConditionsPartitions[i_condition];

        //for each element in the model loop over its connectivities
        for(IO::ConnectivitiesContainerType::value_type::iterator i_node = ConditionsConnectivities[i_condition].begin() ;
                i_node != ConditionsConnectivities[i_condition].end() ; i_node++)
        {
            //get global id. We assume that node ids are began with one
            const int my_gid = *i_node-1;

            //get the partition index for the node i am interested in
            const int node_partition = NodesPartitions[my_gid];

            // adding the partition of the element to its nodes
            if(condition_partition != node_partition) // we will add the node_partition once afterward
                rNodesAllPartitions[my_gid].push_back(condition_partition);
        }
    }

    // adding the nodes partition to their array of partitions and clear the repeated ones
    for(SizeType i_node = 0 ; i_node < number_of_nodes ; i_node++)
    {
        IO::PartitionIndicesContainerType::value_type& node_partitions = rNodesAllPartitions[i_node];
        node_partitions.push_back(NodesPartitions[i_node]);

        std::sort(node_partitions.begin(), node_partitions.end());
        IO::PartitionIndicesContainerType::value_type::iterator new_end=std::unique(node_partitions.begin(), node_partitions.end());
        node_partitions.resize(new_end - node_partitions.begin());
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
