#ifndef KRATOS_METIS_DIVIDE_HETEROGENEOUS_INPUT_PROCESS_H
#define KRATOS_METIS_DIVIDE_HETEROGENEOUS_INPUT_PROCESS_H

// External includes
#include <parmetis.h>

// Project includes
#include "includes/define.h"
#include "includes/io.h"
#include "processes/process.h"
#include "processes/graph_coloring_process.h"
#include "custom_processes/metis_divide_input_to_partitions_process.h"

extern "C" {
extern void METIS_PartGraphKway(int*,  //int* n
                                int*,  //idxtype* xadj
                                int*,  //idxtype* adjcncy
                                int*,  //idxtype* vwgt
                                int*,  //idxtype* adjwgt
                                int*,  //int* wgtflag
                                int*,  //int* numflag
                                int*,  //int* nparts
                                int*,  //int* options
                                int*,  //int* edgecut
                                int*); //indxtype* part
};

namespace Kratos
{
///@addtogroup MetisApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Call Metis to divide an heterogeneous mesh, by partitioning its nodal graph.
class MetisDivideHeterogeneousInputProcess : public MetisDivideInputToPartitionsProcess
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MetisDivideHeterogeneousInputProcess
    KRATOS_CLASS_POINTER_DEFINITION(MetisDivideHeterogeneousInputProcess);

    typedef MetisDivideInputToPartitionsProcess BaseType;

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;
    typedef matrix<int> GraphType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MetisDivideHeterogeneousInputProcess(IO& rIO, SizeType NumberOfPartitions, int Dimension = 3,int Verbosity = 0):
        MetisDivideInputToPartitionsProcess(rIO,NumberOfPartitions,Dimension),
        mVerbosity(Verbosity)
    {
    }

    /// Copy constructor.
    MetisDivideHeterogeneousInputProcess(MetisDivideHeterogeneousInputProcess const& rOther):
        MetisDivideInputToPartitionsProcess(rOther.mrIO,rOther.mNumberOfPartitions,rOther.mDimension),
        mVerbosity(rOther.mVerbosity)
    {
    }

    /// Destructor.
    virtual ~MetisDivideHeterogeneousInputProcess()
    {
    }


    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        this->Execute();
    }


    ///@}
    ///@name Operations
    ///@{

    /// Generate a partition using Metis.
    /** Partitioned input is written as <problem name>_<mpi rank>.mdpa
     */
    virtual void Execute()
    {
        // Read nodal graph from input
        int* NodeIndices = 0;
        int* NodeConnectivities = 0;

        SizeType NumNodes = BaseType::mrIO.ReadNodalGraph(&NodeIndices,&NodeConnectivities);

        std::vector<int> NodePartition;
        PartitionNodes(NumNodes,NodeIndices,NodeConnectivities,NodePartition);

        // Liberate some memory we no longer need
        delete [] NodeIndices;
        delete [] NodeConnectivities;

        // Partition elements
        IO::ConnectivitiesContainerType ElementConnectivities;
        SizeType NumElements =  BaseType::mrIO.ReadElementsConnectivities(ElementConnectivities);

        std::vector<int> ElementPartition;

        PartitionMesh(NodePartition,ElementConnectivities,ElementPartition);

        // Partition conditions
        IO::ConnectivitiesContainerType ConditionConnectivities;
        SizeType NumConditions = BaseType::mrIO.ReadConditionsConnectivities(ConditionConnectivities);

        std::vector<int> ConditionPartition;

        PartitionMesh(NodePartition,ConditionConnectivities,ConditionPartition);

        // Coloring
        GraphType DomainGraph = zero_matrix<int>(mNumberOfPartitions);
        CalculateDomainsGraph(DomainGraph,NumElements,ElementConnectivities,NodePartition,ElementPartition);
        CalculateDomainsGraph(DomainGraph,NumConditions,ConditionConnectivities,NodePartition,ConditionPartition);

        int NumColors;
        GraphType ColoredDomainGraph;
        GraphColoringProcess(mNumberOfPartitions,DomainGraph,ColoredDomainGraph,NumColors).Execute();

        if (mVerbosity > 0)
        {
            KRATOS_WATCH(NumColors);
            KRATOS_WATCH(ColoredDomainGraph);
        }

        // Write partition info into separate input files
        IO::PartitionIndicesContainerType nodes_all_partitions;
        IO::PartitionIndicesContainerType elements_all_partitions;
        IO::PartitionIndicesContainerType conditions_all_partitions;

        // Create lists containing all nodes/elements/conditions known to each partition
        DividingNodes(nodes_all_partitions, ElementConnectivities, ConditionConnectivities, NodePartition, ElementPartition, ConditionPartition);
        DividingElements(elements_all_partitions, ElementPartition);
        DividingConditions(conditions_all_partitions, ConditionPartition);

        if (mVerbosity > 1)
        {
            std::cout << "Final list of nodes known by each partition" << std::endl;
            for(SizeType i = 0 ; i < NumNodes ; i++)
            {
                std::cout << "Node #" << i+1 << "->";
                for(std::vector<std::size_t>::iterator j = nodes_all_partitions[i].begin() ; j != nodes_all_partitions[i].end() ; j++)
                    std::cout << *j << ",";
                std::cout << std::endl;
            }
        }

        IO::PartitionIndicesType io_nodes_partitions(NodePartition.begin(), NodePartition.end());
        IO::PartitionIndicesType io_elements_partitions(ElementPartition.begin(), ElementPartition.end());
        IO::PartitionIndicesType io_conditions_partitions(ConditionPartition.begin(), ConditionPartition.end());

        // Write files
        mrIO.DivideInputToPartitions(mNumberOfPartitions, ColoredDomainGraph,
                                     io_nodes_partitions, io_elements_partitions, io_conditions_partitions,
                                     nodes_all_partitions, elements_all_partitions, conditions_all_partitions);
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "MetisDivideHeterogeneousInputProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MetisDivideHeterogeneousInputProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    int mVerbosity;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /// Call Metis to create a partition of nodes based on the nodal graph.
    /**
     * Metis_PartGraphKway is used to generate a partition for nodes in the model.
     * @param NumNodes Number of nodes in the mesh.
     * @param NodeIndices pointer to C array of size NumNodes+1, xadj array for calls to Metis_PartGraphKway.
     * @param NodeConnectivities pointer to C array of size NodeIndices[NumNodes+1], adjcncy array for calls to Metis_PartGraphKway.
     * @param rNodePartition Metis output as an std::vector, position k is the partition that will contain the k-th node (Node with Id k+1).
     * @return Number of graph edges cut by Metis.
     */
    int PartitionNodes(SizeType NumNodes,
                       int* NodeIndices,
                       int* NodeConnectivities,
                       std::vector<int>& rNodePartition)
    {
        int n = static_cast<int>(NumNodes);
        int wgtflag = 0; // Graph is not weighted
        int numflag = 0; // Nodes are numbered from 0 (C style)
        int nparts = static_cast<int>(BaseType::mNumberOfPartitions);
        int options[5]; // options array
        options[0] = 0; // use default options
        int edgecut;
        rNodePartition.resize(NumNodes);

        METIS_PartGraphKway(&n,NodeIndices,NodeConnectivities,NULL,NULL,&wgtflag,&numflag,&nparts,&options[0],&edgecut,&rNodePartition[0]);

        // Debug: print partition
        PrintDebugData("Node Partition",rNodePartition);

        return edgecut;
    }

    /// Use the nodal partition data to assign elements or conditions to a partition.
    void PartitionMesh(std::vector<int> const& NodePartition,
                       const IO::ConnectivitiesContainerType& rElemConnectivities,
                       std::vector<int>& rElemPartition)
    {
        SizeType NumElements = rElemConnectivities.size();
        std::vector<int> PartitionWeights(BaseType::mNumberOfPartitions,0);

        // initialize ElementPartition
        rElemPartition.resize(NumElements,-1);

        // Elements where all nodes belong to the same partition always go to that partition
        IO::ConnectivitiesContainerType::const_iterator itElem = rElemConnectivities.begin();
        for (std::vector<int>::iterator itPart = rElemPartition.begin(); itPart != rElemPartition.end(); itPart++)
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
        int MaxWeight = 1.03 * NumElements / BaseType::mNumberOfPartitions;
        for (std::vector<int>::iterator itPart = rElemPartition.begin(); itPart != rElemPartition.end(); itPart++)
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

    SizeType FindMax(SizeType NumTerms, const std::vector<int>& rVect)
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

    void PrintDebugData(const std::string& rLabel,
                        const std::vector<int>& rPartitionData)
    {
        if (mVerbosity > 0)
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
                        if (mVerbosity > 1)
                            std::cout << i+1 << ",";
                    }
                }
                std::cout << " contains " << count << " items." << std::endl;
            }
        }
    }


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    MetisDivideHeterogeneousInputProcess& operator=(MetisDivideHeterogeneousInputProcess const& rOther);

    // Copy constructor.
    //MetisDivideHeterogeneousInputProcess(MetisDivideHeterogeneousInputProcess const& rOther);


    ///@}

}; // Class MetisDivideHeterogeneousInputProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  MetisDivideHeterogeneousInputProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MetisDivideHeterogeneousInputProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // addtogroup block

}

#endif // KRATOS_METIS_DIVIDE_HETEROGENEOUS_INPUT_PROCESS_H
