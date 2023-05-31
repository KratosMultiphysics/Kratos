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
    SizeType NumNodes;
    std::vector<idxtype> NodePartition;
    GetNodesPartitions(NodePartition, NumNodes);

    // Partition elements
    IO::ConnectivitiesContainerType ElementConnectivities;
    SizeType NumElements =  mrIO.ReadElementsConnectivities(ElementConnectivities);
    if (NumElements != ElementConnectivities.size())
    {
        std::stringstream Msg;
        Msg << std::endl;
        Msg << "ERROR in MetisDivideHeterogenousInputProcess:" << std::endl;
        Msg << "Read " << NumElements << " elements, but element list has " << ElementConnectivities.size() << " entries." << std::endl;
        Msg << "Elements are most likely not correlatively numbered." << std::endl;

        KRATOS_THROW_ERROR(std::runtime_error,Msg.str(),"");
    }

    std::vector<idxtype> ElementPartition;

    if (mSynchronizeConditions)
        PartitionElementsSynchronous(NodePartition,ElementConnectivities,ElementPartition);
    else
        PartitionMesh(NodePartition,ElementConnectivities,ElementPartition);

    // Partition conditions
    IO::ConnectivitiesContainerType ConditionConnectivities;
    SizeType NumConditions = mrIO.ReadConditionsConnectivities(ConditionConnectivities);
    if (NumConditions != ConditionConnectivities.size())
    {
        std::stringstream Msg;
        Msg << std::endl;
        Msg << "ERROR in MetisDivideHeterogenousInputProcess:" << std::endl;
        Msg << "Read " << NumConditions << " conditions, but condition list has " << ConditionConnectivities.size() << " entries." << std::endl;
        Msg << "Conditions are most likely not correlatively numbered." << std::endl;

        KRATOS_THROW_ERROR(std::runtime_error,Msg.str(),"");
    }

    std::vector<idxtype> ConditionPartition;

    if (mSynchronizeConditions)
        PartitionConditionsSynchronous(NodePartition,ElementPartition,ConditionConnectivities,ElementConnectivities,ConditionPartition);
    else
        PartitionMesh(NodePartition,ConditionConnectivities,ConditionPartition);

    // Detect hanging nodes (nodes that belong to a partition where no local elements have them) and send them to another partition.
    // Hanging nodes should be avoided, as they can cause problems when setting the Dofs
    RedistributeHangingNodes(NodePartition,ElementPartition,ElementConnectivities,ConditionPartition,ConditionConnectivities);

    // Coloring
    GraphType DomainGraph = zero_matrix<int>(mNumberOfPartitions);
    LegacyPartitioningUtilities::CalculateDomainsGraph(DomainGraph,NumElements,ElementConnectivities,NodePartition,ElementPartition);
    LegacyPartitioningUtilities::CalculateDomainsGraph(DomainGraph,NumConditions,ConditionConnectivities,NodePartition,ConditionPartition);

    int NumColors;
    GraphColoringProcess(mNumberOfPartitions,DomainGraph,rPartitioningInfo.Graph,NumColors).Execute();

    if (mVerbosity > 0)
    {
        KRATOS_WATCH(NumColors);
    }

if (mVerbosity > 2)
{
        KRATOS_WATCH(rPartitioningInfo.Graph);
}

    // Write partition info into separate input files
    // Create lists containing all nodes/elements/conditions known to each partition
    LegacyPartitioningUtilities::DividingNodes(rPartitioningInfo.NodesAllPartitions, ElementConnectivities, ConditionConnectivities, NodePartition, ElementPartition, ConditionPartition);
    LegacyPartitioningUtilities::DividingElements(rPartitioningInfo.ElementsAllPartitions, ElementPartition);
    LegacyPartitioningUtilities::DividingConditions(rPartitioningInfo.ConditionsAllPartitions, ConditionPartition);

    if (mVerbosity > 1)
    {
        auto& nodes_all_partitions = rPartitioningInfo.NodesAllPartitions;
        std::cout << "Final list of nodes known by each partition" << std::endl;
        for(SizeType i = 0 ; i < NumNodes ; i++)
        {
            std::cout << "Node #" << i+1 << "->";
            for(std::vector<std::size_t>::iterator j = nodes_all_partitions[i].begin() ; j != nodes_all_partitions[i].end() ; j++)
                std::cout << *j << ",";
            std::cout << std::endl;
        }
    }

    rPartitioningInfo.NodesPartitions.assign(NodePartition.begin(), NodePartition.end());
    rPartitioningInfo.ElementsPartitions.assign(ElementPartition.begin(), ElementPartition.end());
    rPartitioningInfo.ConditionsPartitions.assign(ConditionPartition.begin(), ConditionPartition.end());
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
            std::cout << "Relocating " << HangingNodes.size() << " isolated nodes." << std::endl;
        else
            std::cout << "No isolated nodes found." << std::endl;
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

        if (mVerbosity > 0)
            std::cout << "Sending node " << HangingNodes[n] << " to partition " << Destination << std::endl;

        rNodePartition[ HangingNodes[n]-1 ] = Destination;
    }

    if (mVerbosity > 0 && HangingNodes.size() > 0)
        std::cout << "Relocated " << HangingNodes.size() << " isolated nodes." << std::endl;
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
