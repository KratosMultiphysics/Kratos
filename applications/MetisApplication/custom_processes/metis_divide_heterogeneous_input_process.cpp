//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Jordi Cotela
//                   Carlos Roig
//                   Vicente Mataix Ferrandiz
//

// System includes

// External includes
#include "metis.h"

// Project includes
#include "processes/graph_coloring_process.h"
#include "metis_divide_heterogeneous_input_process.h"
#include "custom_utilities/legacy_partitioning_utilities.h" // TODO remove

namespace Kratos {

void MetisDivideHeterogeneousInputProcess::ExecutePartitioning(PartitioningInfo& rPartitioningInfo)
{
    SizeType number_of_nodes;
    std::vector<idxtype> node_partition;
    GetNodesPartitions(node_partition, number_of_nodes);

    // Partition geometries
    IO::ConnectivitiesContainerType geometry_connectivities;
    SizeType number_of_geometries = mrIO.ReadGeometriesConnectivities(geometry_connectivities);
    if (number_of_geometries != geometry_connectivities.size()) {
        std::stringstream Msg;
        Msg << std::endl;
        Msg << "ERROR in MetisDivideHeterogenousInputProcess:" << std::endl;
        Msg << "Read " << number_of_geometries << " geometries, but geometry list has " << geometry_connectivities.size() << " entries." << std::endl;
        Msg << "Geometries are most likely not correlatively numbered." << std::endl;
        KRATOS_THROW_ERROR(std::runtime_error,Msg.str(),"");
    }

    std::vector<idxtype> geometry_partition;

    if (mSynchronizeConditions)
        PartitionElementsSynchronous(node_partition,geometry_connectivities, geometry_partition);
    else
        PartitionMesh(node_partition, geometry_connectivities, geometry_partition);

    // Partition elements
    IO::ConnectivitiesContainerType element_connectivities;
    SizeType number_of_elements =  mrIO.ReadElementsConnectivities(element_connectivities);
    if (number_of_elements != element_connectivities.size()) {
        std::stringstream Msg;
        Msg << std::endl;
        Msg << "ERROR in MetisDivideHeterogenousInputProcess:" << std::endl;
        Msg << "Read " << number_of_elements << " elements, but element list has " << element_connectivities.size() << " entries." << std::endl;
        Msg << "Elements are most likely not correlatively numbered." << std::endl;
        KRATOS_THROW_ERROR(std::runtime_error,Msg.str(),"");
    }

    std::vector<idxtype> element_partition;

    if (mSynchronizeConditions)
        PartitionElementsSynchronous(node_partition,element_connectivities, element_partition);
    else
        PartitionMesh(node_partition, element_connectivities, element_partition);

    // Partition conditions
    IO::ConnectivitiesContainerType condition_connectivities;
    SizeType number_of_conditions = mrIO.ReadConditionsConnectivities(condition_connectivities);
    if (number_of_conditions != condition_connectivities.size()) {
        std::stringstream Msg;
        Msg << std::endl;
        Msg << "ERROR in MetisDivideHeterogenousInputProcess:" << std::endl;
        Msg << "Read " << number_of_conditions << " conditions, but condition list has " << condition_connectivities.size() << " entries." << std::endl;
        Msg << "Conditions are most likely not correlatively numbered." << std::endl;
        KRATOS_THROW_ERROR(std::runtime_error,Msg.str(),"");
    }
    std::vector<idxtype> condition_partition;
    if (mSynchronizeConditions)
        PartitionConditionsSynchronous(node_partition,element_partition,condition_connectivities,element_connectivities,condition_partition);
    else
        PartitionMesh(node_partition,condition_connectivities,condition_partition);

    // Partition master-slave constraints // TODO
    IO::ConnectivitiesContainerType master_slave_constraints_connectivities;
    SizeType number_of_master_slave_constraints = mrIO.ReadMasterSlaveConstraintsConnectivities(master_slave_constraints_connectivities);
    if (number_of_master_slave_constraints != master_slave_constraints_connectivities.size()) {
        std::stringstream Msg;
        Msg << std::endl;
        Msg << "ERROR in MetisDivideHeterogenousInputProcess:" << std::endl;
        Msg << "Read " << number_of_master_slave_constraints << " master-slave constraints, but constraint list has " << master_slave_constraints_connectivities.size() << " entries." << std::endl;
        Msg << "Master-slave constraints are most likely not correlatively numbered." << std::endl;
        KRATOS_THROW_ERROR(std::runtime_error,Msg.str(),"");
    }
    std::vector<idxtype> master_slave_constraints_partition;
    if (mSynchronizeMasterSlaveConstraints)
        PartitionMasterSlaveConstraintsSynchronous(node_partition,master_slave_constraints_connectivities,master_slave_constraints_partition);
    else
        PartitionMesh(node_partition,master_slave_constraints_connectivities,master_slave_constraints_partition);

    // Detect hanging nodes (nodes that belong to a partition where no local elements have them) and send them to another partition.
    // Hanging nodes should be avoided, as they can cause problems when setting the Dofs
    RedistributeHangingNodes(node_partition,element_partition,element_connectivities,condition_partition,condition_connectivities);

    // Coloring
    GraphType DomainGraph = zero_matrix<int>(mNumberOfPartitions);
    LegacyPartitioningUtilities::CalculateDomainsGraph(DomainGraph,number_of_elements,element_connectivities,node_partition,element_partition);
    LegacyPartitioningUtilities::CalculateDomainsGraph(DomainGraph,number_of_conditions,condition_connectivities,node_partition,condition_partition);

    int number_of_colors;
    GraphColoringProcess(mNumberOfPartitions,DomainGraph,rPartitioningInfo.Graph,number_of_colors).Execute();

    KRATOS_INFO_IF("MetisDivideHeterogeneousInputProcess", mVerbosity > 0) << "Number of colors: " << number_of_colors << std::endl;

    KRATOS_INFO_IF("MetisDivideHeterogeneousInputProcess", mVerbosity > 2) << "Graph: " << rPartitioningInfo.Graph << std::endl;

    // Write partition info into separate input files
    // Create lists containing all nodes/elements/conditions known to each partition
    LegacyPartitioningUtilities::DividingNodes(rPartitioningInfo.NodesAllPartitions, element_connectivities, condition_connectivities, node_partition, element_partition, condition_partition);
    LegacyPartitioningUtilities::DividingGeometries(rPartitioningInfo.GeometriesAllPartitions, geometry_partition);
    LegacyPartitioningUtilities::DividingElements(rPartitioningInfo.ElementsAllPartitions, element_partition);
    LegacyPartitioningUtilities::DividingConditions(rPartitioningInfo.ConditionsAllPartitions, condition_partition);
    LegacyPartitioningUtilities::DividingMasterSlaveConstraints(rPartitioningInfo.MasterSlaveConstraintsAllPartitions, master_slave_constraints_partition);

    if (mVerbosity > 1) {
        auto& r_nodes_all_partitions = rPartitioningInfo.NodesAllPartitions;
        std::cout << "Final list of nodes known by each partition" << std::endl;
        for(SizeType i = 0 ; i < number_of_nodes ; i++) {
            std::cout << "Node #" << i+1 << "->";
            for(std::vector<std::size_t>::iterator j = r_nodes_all_partitions[i].begin() ; j != r_nodes_all_partitions[i].end() ; j++)
                std::cout << *j << ",";
            std::cout << std::endl;
        }
    }

    rPartitioningInfo.NodesPartitions.assign(node_partition.begin(), node_partition.end());
    rPartitioningInfo.GeometriesPartitions.assign(geometry_partition.begin(), geometry_partition.end());
    rPartitioningInfo.ElementsPartitions.assign(element_partition.begin(), element_partition.end());
    rPartitioningInfo.ConditionsPartitions.assign(condition_partition.begin(), condition_partition.end());
    rPartitioningInfo.MasterSlaveConstraints.assign(master_slave_constraints_partition.begin(), master_slave_constraints_partition.end());
}

void MetisDivideHeterogeneousInputProcess::Execute()
{
    PartitioningInfo part_info;

    ExecutePartitioning(part_info);

    // Write files
    mrIO.DivideInputToPartitions(
        mNumberOfPartitions,
        part_info);
}


void MetisDivideHeterogeneousInputProcess::GetNodesPartitions(std::vector<idxtype> &rNodePartition, SizeType &rNumNodes)
{
    // Read nodal graph from input
    IO::ConnectivitiesContainerType kratos_format_node_connectivities;

    rNumNodes = mrIO.ReadNodalGraph(kratos_format_node_connectivities);

    SizeType num_nodes_in_mesh = mrIO.ReadNodesNumber();
    if (rNumNodes != num_nodes_in_mesh)
        KRATOS_ERROR << "Invalid mesh: number of connected nodes = " << rNumNodes
                        << ", number of mesh nodes = " << num_nodes_in_mesh << "."
                        << std::endl;

    // Write connectivity data in CSR format
    idxtype* node_indices = 0;
    idxtype* node_connectivities = 0;

    LegacyPartitioningUtilities::ConvertKratosToCSRFormat(kratos_format_node_connectivities, &node_indices, &node_connectivities);

    PartitionNodes(rNumNodes, node_indices, node_connectivities, rNodePartition);

    // Free some memory we no longer need
    delete [] node_indices;
    delete [] node_connectivities;
}

int MetisDivideHeterogeneousInputProcess::PartitionNodes(SizeType NumNodes,
                    idxtype* NodeIndices,
                    idxtype* NodeConnectivities,
                    std::vector<idxtype>& rNodePartition)
{
    mNumNodes = NumNodes;
    idxtype n = static_cast<idxtype>(NumNodes);

    idxtype nparts = static_cast<idxtype>(mNumberOfPartitions);
    idxtype edgecut;
    rNodePartition.resize(NumNodes);

    idxtype ncon = 1; //The number of balancing constraints. It should be at least 1.

    idxtype options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);

    int metis_return = METIS_PartGraphKway(&n,&ncon,NodeIndices,NodeConnectivities,NULL,NULL,NULL,&nparts,NULL,NULL,&options[0],&edgecut,&rNodePartition[0]);

    if(metis_return != METIS_OK)
        std::cout << "metis returns the following error code :" << metis_return << std::endl;
    /*         int METIS PartGraphKway(
        * idx_t *nvtxs,
        * idx_t *ncon,
        * idx_t *xadj,
        * idx_t *adjncy,
        * idx_t *vwgt,   NULL
        * idx_t *vsize,  NULL
        * idx_t *adjwgt,  NULL
        * idx_t *nparts,
        * real_t *tpwgts, NULL
        * real_t ubvec, NULL
        * idx_t *options,
        * idx_t *objval, idx t *part)
        */

    // Debug: print partition
    PrintDebugData("Node Partition",rNodePartition);

    return edgecut;
}

void MetisDivideHeterogeneousInputProcess::PartitionMesh(std::vector<idxtype> const& NodePartition,
                    const IO::ConnectivitiesContainerType& rElemConnectivities,
                    std::vector<idxtype>& rElemPartition)
{
    SizeType NumElements = rElemConnectivities.size();
    std::vector<int> PartitionWeights(mNumberOfPartitions,0);

    // initialize ElementPartition
    rElemPartition.resize(NumElements,-1);

    // Elements where all nodes belong to the same partition always go to that partition
    IO::ConnectivitiesContainerType::const_iterator itElem = rElemConnectivities.begin();
    for (std::vector<idxtype>::iterator itPart = rElemPartition.begin(); itPart != rElemPartition.end(); itPart++)
    {
        int MyPartition = NodePartition[ (*itElem)[0] - 1 ]; // Node Ids start from 1
        SizeType NeighbourNodes = 1; // Nodes in the same partition
        for (std::vector<SizeType>::const_iterator itNode = itElem->begin()+1; itNode != itElem->end(); ++itNode)
        {
            if ( NodePartition[ *itNode - 1 ] == MyPartition )
                ++NeighbourNodes;
            else
                break;
        }

        if ( NeighbourNodes == itElem->size() )
        {
            *itPart = MyPartition;
            PartitionWeights[MyPartition]++;
        }

        // Advance to next element in connectivities array
        itElem++;
    }

    // Now distribute boundary elements
    itElem = rElemConnectivities.begin();
    int MaxWeight = 1.03 * NumElements / mNumberOfPartitions;
    for (std::vector<idxtype>::iterator itPart = rElemPartition.begin(); itPart != rElemPartition.end(); itPart++)
    {
        if (*itPart == -1) // If element is still unassigned
        {
            SizeType FoundNeighbours = 0;
            SizeType NodesInElem = itElem->size();
            std::vector<int> NeighbourPartitions(NodesInElem,-1);
            std::vector<int> NeighbourWeights(NodesInElem,0);

            for (std::vector<SizeType>::const_iterator itNode = itElem->begin(); itNode != itElem->end(); ++itNode)
            {
                // Check if the node's partition was already found in this element
                int MyPartition = NodePartition[ *itNode - 1 ]; // This node's partition
                SizeType i=0;
                for (i = 0; i < FoundNeighbours; i++)
                {
                    if (MyPartition == NeighbourPartitions[i])
                    {
                        NeighbourWeights[i]++;
                        break;
                    }
                }

                // If this is the first node in this partition, add the partition to the candidate partition list
                if (i == FoundNeighbours)
                {
                    NeighbourWeights[i] = 1;
                    NeighbourPartitions[i] = MyPartition;
                    FoundNeighbours++;
                }
            }

            // Determine the partition that owns the most nodes, and try to assign the element to that partition
            // (if the partition is full, try the other partitions)
            int MajorityPartition = NeighbourPartitions[ FindMax(FoundNeighbours,NeighbourWeights) ];
            if(PartitionWeights[MajorityPartition] < MaxWeight)
            {
                *itPart = MajorityPartition;
                PartitionWeights[MajorityPartition]++;
            }
            else
            {
                SizeType i=0;
                for(i = 0; i < FoundNeighbours; i++)
                {
                    int DestPartition = NeighbourPartitions[i];
                    if(PartitionWeights[ DestPartition ] < MaxWeight)
                    {
                        *itPart = DestPartition;
                        PartitionWeights[ DestPartition ]++;
                        break;
                    }
                }

                // If all partitions are full, assign the element to the partition that owns the most nodes
                if (i == FoundNeighbours)
                {
                    *itPart = MajorityPartition;
                    PartitionWeights[MajorityPartition]++;
                }
            }
        }

        // Advance to next element in connectivities array
        itElem++;
    }

    PrintDebugData("Mesh Partition",rElemPartition);

}

void MetisDivideHeterogeneousInputProcess::PartitionGeometriesSynchronous(
    const std::vector<idxtype>& rNodePartition,
    const IO::ConnectivitiesContainerType& rGeomConnectivities,
    std::vector<idxtype>& rGeomPartition
    )
{
    const SizeType number_of_geometries = rGeomConnectivities.size();

    // initialize Geometry Partition
    rGeomPartition.resize(number_of_geometries, -1);

    // Geometries where all nodes belong to the same partition always go to that partition
    auto it_geom = rGeomConnectivities.begin();
    for (auto it_part = rGeomPartition.begin(); it_part != rGeomPartition.end(); it_part++) {
        const int my_partition = rNodePartition[(*it_geom)[0] - 1]; // Node Ids start from 1
        SizeType neighbour_nodes = 1; // Nodes in the same partition
        for (auto it_node = it_geom->begin() + 1; it_node != it_geom->end(); ++it_node) {
            if (rNodePartition[*it_node - 1] == my_partition)
                ++neighbour_nodes;
            else
                break;
        }

        if (neighbour_nodes == it_geom->size()) {
            *it_part = my_partition;
        }

        // Advance to next geometry in the connectivities array
        it_geom++;
    }

    // Now distribute boundary geometries
    it_geom = rGeomConnectivities.begin();
    for (auto it_part = rGeomPartition.begin(); it_part != rGeomPartition.end(); it_part++) {
        if (*it_part == -1) { // If geometry is still unassigned
            SizeType found_neighbours = 0;
            SizeType nodes_in_geom = it_geom->size();
            std::vector<int> neighbour_partitions(nodes_in_geom, -1);
            std::vector<int> neighbour_weights(nodes_in_geom, 0);

            for (auto it_node = it_geom->begin(); it_node != it_geom->end(); ++it_node) {
                // Check if the node's partition was already found in this geometry
                const int my_partition = rNodePartition[*it_node - 1]; // This node's partition
                SizeType i = 0;
                for (i = 0; i < found_neighbours; i++) {
                    if (my_partition == neighbour_partitions[i]) {
                        neighbour_weights[i]++;
                        break;
                    }
                }

                // If this is the first node in this partition, add the partition to the candidate list
                if (i == found_neighbours) {
                    neighbour_weights[i] = 1;
                    neighbour_partitions[i] = my_partition;
                    found_neighbours++;
                }
            }

            // Determine the partition that owns the most nodes
            const int majority_partition = neighbour_partitions[FindMax(found_neighbours, neighbour_weights)];
            *it_part = majority_partition;
        }

        // Advance to next geometry in the connectivities array
        it_geom++;
    }

    PrintDebugData("Geometry Partition", rGeomPartition);
}

void MetisDivideHeterogeneousInputProcess::PartitionElementsSynchronous(std::vector<idxtype> const& NodePartition,
                    const IO::ConnectivitiesContainerType& rElemConnectivities,
                    std::vector<idxtype>& rElemPartition)
{
    SizeType NumElements = rElemConnectivities.size();

    // initialize ElementPartition
    mNodeConnectivities = std::vector<std::unordered_set<std::size_t>>(mNumNodes,std::unordered_set<std::size_t>());
    rElemPartition.resize(NumElements,-1);

    // Fill the node Connectivities
    for(std::size_t i = 0; i < NumElements; i++) {
        for (std::vector<SizeType>::const_iterator itNode = rElemConnectivities[i].begin(); itNode != rElemConnectivities[i].end(); ++itNode) {
            mNodeConnectivities[*itNode-1].insert(i);
        }
    }

    // Elements where all nodes belong to the same partition always go to that partition
    IO::ConnectivitiesContainerType::const_iterator itElem = rElemConnectivities.begin();
    for (std::vector<idxtype>::iterator itPart = rElemPartition.begin(); itPart != rElemPartition.end(); itPart++)
    {
        int MyPartition = NodePartition[ (*itElem)[0] - 1 ]; // Node Ids start from 1
        SizeType NeighbourNodes = 1; // Nodes in the same partition
        for (std::vector<SizeType>::const_iterator itNode = itElem->begin()+1; itNode != itElem->end(); ++itNode)
        {
            if ( NodePartition[ *itNode - 1 ] == MyPartition )
                ++NeighbourNodes;
            else
                break;
        }

        if ( NeighbourNodes == itElem->size() )
        {
            *itPart = MyPartition;
        }

        // Advance to next element in connectivities array
        itElem++;
    }

    // Now distribute boundary elements
    itElem = rElemConnectivities.begin();
    //int MaxWeight = 1.03 * NumElements / mNumberOfPartitions;
    for (std::vector<idxtype>::iterator itPart = rElemPartition.begin(); itPart != rElemPartition.end(); itPart++)
    {
        if (*itPart == -1) // If element is still unassigned
        {
            SizeType FoundNeighbours = 0;
            SizeType NodesInElem = itElem->size();
            std::vector<int> NeighbourPartitions(NodesInElem,-1);
            std::vector<int> NeighbourWeights(NodesInElem,0);

            for (std::vector<SizeType>::const_iterator itNode = itElem->begin(); itNode != itElem->end(); ++itNode)
            {
                // Check if the node's partition was already found in this element
                int MyPartition = NodePartition[ *itNode - 1 ]; // This node's partition
                SizeType i=0;
                for (i = 0; i < FoundNeighbours; i++)
                {
                    if (MyPartition == NeighbourPartitions[i])
                    {
                        NeighbourWeights[i]++;
                        break;
                    }
                }

                // If this is the first node in this partition, add the partition to the candidate partition list
                if (i == FoundNeighbours)
                {
                    NeighbourWeights[i] = 1;
                    NeighbourPartitions[i] = MyPartition;
                    FoundNeighbours++;
                }
            }

            // Determine the partition that owns the most nodes, and try to assign the element to that partition
            int MajorityPartition = NeighbourPartitions[ FindMax(FoundNeighbours,NeighbourWeights) ];
            {
                *itPart = MajorityPartition;
            }
        }

        // Advance to next element in connectivities array
        itElem++;
    }

    PrintDebugData("Element Partition",rElemPartition);

}

void MetisDivideHeterogeneousInputProcess::PartitionConditionsSynchronous(const std::vector<idxtype>& rNodePartition,
                const std::vector<idxtype>& rElemPartition,
                const IO::ConnectivitiesContainerType& rCondConnectivities,
                const IO::ConnectivitiesContainerType& rElemConnectivities,
                std::vector<idxtype>& rCondPartition)
{
    SizeType NumElements = rElemConnectivities.size();
    SizeType NumConditions = rCondConnectivities.size();

    // initialize CondPartition
    rCondPartition.resize(NumConditions,-1);

    // make sorted element connectivities array
    IO::ConnectivitiesContainerType ElementsSorted(rElemConnectivities);
    for (SizeType i=0; i<NumElements; i++)
        std::sort(ElementsSorted[i].begin(), ElementsSorted[i].end());

    IO::ConnectivitiesContainerType::const_iterator itCond = rCondConnectivities.begin();
    //int MaxWeight = 1.03 * NumConditions / mNumberOfPartitions;
    for (std::vector<idxtype>::iterator itPart = rCondPartition.begin(); itPart != rCondPartition.end(); itPart++)
    {
        SizeType FoundNeighbours = 0;
        SizeType NodesInCond = itCond->size();
        std::vector<int> NeighbourPartitions(NodesInCond,-1);
        std::vector<int> NeighbourWeights(NodesInCond,0);

        for (std::vector<SizeType>::const_iterator itNode = itCond->begin(); itNode != itCond->end(); ++itNode)
        {
            // Check if the node's partition was already found in this condition
            int MyPartition = rNodePartition[ *itNode - 1 ]; // This node's partition
            SizeType i=0;
            for (i = 0; i < FoundNeighbours; i++)
            {
                if (MyPartition == NeighbourPartitions[i])
                {
                    NeighbourWeights[i]++;
                    break;
                }
            }

            // If this is the first node in this partition, add the partition to the candidate partition list
            if (i == FoundNeighbours)
            {
                NeighbourWeights[i] = 1;
                NeighbourPartitions[i] = MyPartition;
                FoundNeighbours++;
            }
        }
        // Determine the partition that owns the most nodes, and try to assign the condition to that partition
        int MajorityPartition = NeighbourPartitions[ FindMax(FoundNeighbours,NeighbourWeights) ];
        {
            *itPart = MajorityPartition;
        }

        // ensure conditions sharing nodes with an element have same partition as the element
        IO::ConnectivitiesContainerType::value_type tmp(*itCond);
        std::sort(tmp.begin(), tmp.end());

        for (std::vector<SizeType>::const_iterator itNode = itCond->begin(); itNode != itCond->end(); ++itNode) {
            for(auto shared_element : mNodeConnectivities[*itNode - 1]) {
                // Would it be faster to sort the element here as well?
                // Should be as far as "numConditions * conPerNode >> numElements", but not otherwise
                if ( std::includes(ElementsSorted[shared_element].begin(), ElementsSorted[shared_element].end(), tmp.begin(), tmp.end()) ) {
                *itPart = rElemPartition[shared_element];
                break;
                }
            }
        }

        // Advance to next condition in connectivities array
        ++itCond;
    }

    PrintDebugData("Condition Partition",rCondPartition);

}

void MetisDivideHeterogeneousInputProcess::PartitionMasterSlaveConstraintsSynchronous(std::vector<idxtype> const& NodePartition,
                    const IO::ConnectivitiesContainerType& rMasterSlaveConstraintConnectivities,
                    std::vector<idxtype>& rMasterSlaveConstraintPartition)
{
    SizeType NumMasterSlaveConstraints = rMasterSlaveConstraintConnectivities.size();

    // initialize MasterSlaveConstraintPartition
    mNodeConnectivities = std::vector<std::unordered_set<std::size_t>>(mNumNodes,std::unordered_set<std::size_t>());
    rMasterSlaveConstraintPartition.resize(NumMasterSlaveConstraints,-1);

    // Fill the node Connectivities
    for(std::size_t i = 0; i < NumMasterSlaveConstraints; i++) {
        for (auto it_node = rMasterSlaveConstraintConnectivities[i].begin(); it_node != rMasterSlaveConstraintConnectivities[i].end(); ++it_node) {
            mNodeConnectivities[*it_node-1].insert(i);
        }
    }

    // MasterSlaveConstraints where all nodes belong to the same partition always go to that partition
    IO::ConnectivitiesContainerType::const_iterator it_master_slave_constraint = rMasterSlaveConstraintConnectivities.begin();
    for (auto it_part = rMasterSlaveConstraintPartition.begin(); it_part != rMasterSlaveConstraintPartition.end(); it_part++) {
        const int my_partition = NodePartition[ (*it_master_slave_constraint)[0] - 1 ]; // Node Ids start from 1
        SizeType neighbour_nodes = 1; // Nodes in the same partition
        for (auto it_node = it_master_slave_constraint->begin()+1; it_node != it_master_slave_constraint->end(); ++it_node) {
            if ( NodePartition[ *it_node - 1 ] == my_partition )
                ++neighbour_nodes;
            else
                break;
        }

        if ( neighbour_nodes == it_master_slave_constraint->size() ) {
            *it_part = my_partition;
        }

        // Advance to next constraint in connectivities array
        it_master_slave_constraint++;
    }

    // Now distribute boundary constraints
    it_master_slave_constraint = rMasterSlaveConstraintConnectivities.begin();
    for (auto it_part = rMasterSlaveConstraintPartition.begin(); it_part != rMasterSlaveConstraintPartition.end(); it_part++) {
        if (*it_part == -1) // If constraint is still unassigned {
            SizeType found_neighbours = 0;
            SizeType nodes_in_master_slave_constraint = it_master_slave_constraint->size();
            std::vector<int> neighbour_partitions(nodes_in_master_slave_constraint,-1);
            std::vector<int> neighbour_weights(nodes_in_master_slave_constraint,0);

            for (auto it_node = it_master_slave_constraint->begin(); it_node != it_master_slave_constraint->end(); ++it_node) {
                // Check if the node's partition was already found in this constraint
                const int my_partition = NodePartition[ *it_node - 1 ]; // This node's partition
                SizeType i=0;
                for (i = 0; i < found_neighbours; i++) {
                    if (my_partition == neighbour_partitions[i]) {
                        neighbour_weights[i]++;
                        break;
                    }
                }

                // If this is the first node in this partition, add the partition to the candidate partition list
                if (i == found_neighbours) {
                    neighbour_weights[i] = 1;
                    neighbour_partitions[i] = my_partition;
                    found_neighbours++;
                }
            }

            // Determine the partition that owns the most nodes, and try to assign the constraint to that partition
            const int majority_partitition = neighbour_partitions[ FindMax(found_neighbours,neighbour_weights) ];
            {
                *it_part = majority_partitition;
            }
        }

        // Advance to next constraint in connectivities array
        it_master_slave_constraint++;
    }

    PrintDebugData("MasterSlaveConstraint Partition",rMasterSlaveConstraintPartition);
}

void MetisDivideHeterogeneousInputProcess::RedistributeHangingNodes(
        std::vector<idxtype>& rNodePartition,
        std::vector<idxtype> const& rElementPartition,
        const IO::ConnectivitiesContainerType& rElementConnectivities,
        std::vector<idxtype> const& rConditionPartition,
        const IO::ConnectivitiesContainerType& rConditionConnectivities)
{
    std::vector<int> NodeUseCounts(rNodePartition.size(),0);

    // Count number of times a node is used locally
    unsigned int ElemIndex = 0;
    for (IO::ConnectivitiesContainerType::const_iterator iElem = rElementConnectivities.begin(); iElem != rElementConnectivities.end(); iElem++)
    {
        for (std::vector<std::size_t>::const_iterator iNode = iElem->begin(); iNode != iElem->end(); iNode++)
            if ( rNodePartition[ *iNode-1] == rElementPartition[ ElemIndex ] )
                NodeUseCounts[ *iNode-1 ]++;
        ElemIndex++;
    }

    unsigned int CondIndex = 0;
    for (IO::ConnectivitiesContainerType::const_iterator iCond = rConditionConnectivities.begin(); iCond != rConditionConnectivities.end(); iCond++)
    {
        for (std::vector<std::size_t>::const_iterator iNode = iCond->begin(); iNode != iCond->end(); iNode++)
            if ( rNodePartition[ *iNode-1] == rConditionPartition[ CondIndex ] )
                NodeUseCounts[ *iNode-1 ]++;
        CondIndex++;
    }

    std::vector<std::size_t> HangingNodes;
    for (unsigned int i = 0; i < NodeUseCounts.size(); i++)
        if( NodeUseCounts[i] == 0 )
            HangingNodes.push_back( i+1 );

    if (mVerbosity > 0)
    {
        if (HangingNodes.size() > 0)
            KRATOS_INFO("MetisDivideHeterogeneousInputProcess") << "Relocating " << HangingNodes.size() << " isolated nodes." << std::endl;
        else
            KRATOS_INFO("MetisDivideHeterogeneousInputProcess") << "No isolated nodes found." << std::endl;
    }

    // Find a new home for hanging nodes
    for (unsigned int n = 0; n < HangingNodes.size(); n++)
    {
        std::vector<int> LocalUseCount(mNumberOfPartitions,0);
        unsigned int ElemIndex = 0;

        for (IO::ConnectivitiesContainerType::const_iterator iElem = rElementConnectivities.begin(); iElem != rElementConnectivities.end(); iElem++)
        {
            for (std::vector<std::size_t>::const_iterator iNode = iElem->begin(); iNode != iElem->end(); iNode++)
                if ( HangingNodes[n] == *iNode )
                    LocalUseCount[ rElementPartition[ElemIndex] ]++;
            ElemIndex++;
        }

        unsigned int CondIndex = 0;
        for (IO::ConnectivitiesContainerType::const_iterator iCond = rConditionConnectivities.begin(); iCond != rConditionConnectivities.end(); iCond++)
        {
            for (std::vector<std::size_t>::const_iterator iNode = iCond->begin(); iNode != iCond->end(); iNode++)
                if ( HangingNodes[n] == *iNode )
                    LocalUseCount[ rConditionPartition[CondIndex] ]++;
            CondIndex++;
        }

        SizeType Destination = FindMax(mNumberOfPartitions,LocalUseCount);

        KRATOS_INFO_IF("MetisDivideHeterogeneousInputProcess", mVerbosity > 0) << "Sending node " << HangingNodes[n] << " to partition " << Destination << std::endl;

        rNodePartition[ HangingNodes[n]-1 ] = Destination;
    }

    KRATOS_INFO_IF("MetisDivideHeterogeneousInputProcess", mVerbosity > 0 && HangingNodes.size() > 0)  << "Relocated " << HangingNodes.size() << " isolated nodes." << std::endl;
}

void MetisDivideHeterogeneousInputProcess::RedistributeHangingNodes(
    std::vector<idxtype>& rNodePartition,
    std::vector<idxtype> const& rGeometryPartition,
    const IO::ConnectivitiesContainerType& rGeometryConnectivities,
    std::vector<idxtype> const& rElementPartition,
    const IO::ConnectivitiesContainerType& rElementConnectivities,
    std::vector<idxtype> const& rConditionPartition,
    const IO::ConnectivitiesContainerType& rConditionConnectivities,
    std::vector<idxtype> const& rMasterSlaveConstraintPartition,
    const IO::ConnectivitiesContainerType& rMasterSlaveConstraintConnectivities
    )
{
    // Initialize node use counts
    std::vector<int> node_use_counts(rNodePartition.size(),0);

    // Count number of times a node is used locally in geometries
    std::size_t geometry_index = 0;
    for (auto it_geom = rGeometryConnectivities.begin(); it_geom != rGeometryConnectivities.end(); it_geom++) {
        for (auto it_node = it_geom->begin(); it_node != it_geom->end(); it_node++) {
            if ( rNodePartition[ *it_node-1] == rGeometryPartition[ geometry_index ] ) {
                node_use_counts[ *it_node-1 ]++;
            }
        }
        geometry_index++;
    }

    // Count number of times a node is used locally in elements
    std::size_t element_index = 0;
    for (auto it_elem = rElementConnectivities.begin(); it_elem != rElementConnectivities.end(); it_elem++) {
        for (auto it_node = it_elem->begin(); it_node != it_elem->end(); it_node++) {
            if ( rNodePartition[ *it_node-1] == rElementPartition[ element_index ] ) {
                node_use_counts[ *it_node-1 ]++;
            }
        }
        element_index++;
    }

    // Count number of times a node is used locally in conditions
    std::size_t condition_index = 0;
    for (auto it_cond = rConditionConnectivities.begin(); it_cond != rConditionConnectivities.end(); it_cond++) {
        for (auto it_node = it_cond->begin(); it_node != it_cond->end(); it_node++) {
            if ( rNodePartition[ *it_node-1] == rConditionPartition[ condition_index ] ) {
                node_use_counts[ *it_node-1 ]++;
            }
        }
        condition_index++;
    }

    // Count number of times a node is used locally in master-slave constraints
    std::size_t constraint_index = 0;
    for (auto it_const = rMasterSlaveConstraintConnectivities.begin(); it_const != rMasterSlaveConstraintConnectivities.end(); it_const++) {
        for (auto it_node = it_const->begin(); it_node != it_const->end(); it_node++) {
            if ( rNodePartition[*it_node-1] == rMasterSlaveConstraintPartition[constraint_index] ) {
                node_use_counts[*it_node-1]++;
            }
        }
        constraint_index++;
    }

    std::vector<std::size_t> hanging_nodes;
    for (std::size_t i = 0; i < node_use_counts.size(); i++) {
        if( node_use_counts[i] == 0 ) {
            hanging_nodes.push_back( i+1 );
        }
    }

    if (mVerbosity > 0) {
        if (hanging_nodes.size() > 0) {
            KRATOS_INFO("MetisDivideHeterogeneousInputProcess") << "Relocating " << hanging_nodes.size() << " isolated nodes." << std::endl;
        } else {
            KRATOS_INFO("MetisDivideHeterogeneousInputProcess") << "No isolated nodes found." << std::endl;
        }
    }

    // Find a new home for hanging nodes
    for (std::size_t n = 0; n < hanging_nodes.size(); n++) {
        std::vector<int> local_use_count(mNumberOfPartitions,0);

        // Count number of times a node is used locally in geometries
        std::size_t geometry_index = 0;
        for (auto it_geom = rGeometryConnectivities.begin(); it_geom != rGeometryConnectivities.end(); it_geom++) {
            for (auto it_node = it_geom->begin(); it_node != it_geom->end(); it_node++) {
                if ( hanging_nodes[n] == *it_node ) {
                    local_use_count[rGeometryPartition[geometry_index]]++;
                }
            }
            geometry_index++;
        }

        // Count number of times a node is used locally in elements
        std::size_t element_index = 0;
        for (auto it_elem = rElementConnectivities.begin(); it_elem != rElementConnectivities.end(); it_elem++) {
            for (auto it_node = it_elem->begin(); it_node != it_elem->end(); it_node++) {
                if ( hanging_nodes[n] == *it_node ) {
                    local_use_count[rElementPartition[element_index]]++;
                }
            }
            element_index++;
        }

        // Count number of times a node is used locally in conditions
        std::size_t condition_index = 0;
        for (auto it_cond = rConditionConnectivities.begin(); it_cond != rConditionConnectivities.end(); it_cond++) {
            for (auto it_node = it_cond->begin(); it_node != it_cond->end(); it_node++) {
                if ( hanging_nodes[n] == *it_node ) {
                    local_use_count[rConditionPartition[condition_index]]++;
                }
            }
            condition_index++;
        }

        // Count number of times a node is used locally in master-slave constraints
        std::size_t constraint_index = 0;
        for (auto it_const = rMasterSlaveConstraintConnectivities.begin(); it_const != rMasterSlaveConstraintConnectivities.end(); it_const++) {
            for (auto it_node = it_const->begin(); it_node != it_const->end(); it_node++) {
                if ( hanging_nodes[n] == *it_node ) {
                    local_use_count[rMasterSlaveConstraintPartition[constraint_index]]++;
                }
            }
            constraint_index++;
        }

        const SizeType destination = FindMax(mNumberOfPartitions,local_use_count);

        KRATOS_INFO_IF("MetisDivideHeterogeneousInputProcess", mVerbosity > 0) << "Sending node " << hanging_nodes[n] << " to partition " << destination << std::endl;

        rNodePartition[hanging_nodes[n]-1] = destination;
    }

    KRATOS_INFO_IF("MetisDivideHeterogeneousInputProcess", mVerbosity > 0 && hanging_nodes.size() > 0)  << "Relocated " << hanging_nodes.size() << " isolated nodes." << std::endl;
}

MetisDivideHeterogeneousInputProcess::SizeType MetisDivideHeterogeneousInputProcess::FindMax(SizeType NumTerms, const std::vector<int>& rVect)
{
    int max = rVect[0];
    SizeType imax = 0;
    for (SizeType i = 1; i < NumTerms; i++)
        if( rVect[i] > max)
        {
            max = rVect[i];
            imax = i;
        }
    return imax;
}

void MetisDivideHeterogeneousInputProcess::PrintDebugData(const std::string& rLabel,
                    const std::vector<idxtype>& rPartitionData)
{
    if (mVerbosity > 1)
    {
        std::cout << rLabel << std::endl;
        for (int p = 0; p < static_cast<int>(mNumberOfPartitions); p++)
        {
            int count = 0;
            std::cout << "Partition " << p << ": ";
            for (SizeType i = 0; i < rPartitionData.size(); i++)
            {
                if(rPartitionData[i] == p)
                {
                    count++;
                    if (mVerbosity > 2)
                        std::cout << i+1 << ",";
                }
            }
            std::cout << count << " objects." << std::endl;
        }
    }
}

} // namespace Kratos
