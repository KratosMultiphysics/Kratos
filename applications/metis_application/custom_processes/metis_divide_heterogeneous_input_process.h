#ifndef KRATOS_METIS_DIVIDE_HETEROGENEOUS_INPUT_PROCESS_H
#define KRATOS_METIS_DIVIDE_HETEROGENEOUS_INPUT_PROCESS_H

#ifdef KRATOS_USE_METIS_5
  #include "metis.h"
#else
  #include <parmetis.h>
#endif

// Project includes
#include "includes/define.h"
#include "includes/io.h"
#include "processes/process.h"
#include "processes/graph_coloring_process.h"
#include "custom_processes/metis_divide_input_to_partitions_process.h"

#ifndef KRATOS_USE_METIS_5
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
  }
#endif

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

    #ifdef KRATOS_USE_METIS_5
      typedef idx_t idxtype;
    #else
      typedef int idxtype;
    #endif

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
    MetisDivideHeterogeneousInputProcess(IO& rIO, SizeType NumberOfPartitions, int Dimension = 3, int Verbosity = 0, bool SynchronizeConditions = false):
        MetisDivideInputToPartitionsProcess(rIO,NumberOfPartitions,Dimension),
        mSynchronizeConditions(SynchronizeConditions),
    mVerbosity(Verbosity)
    {
    }

    /// Copy constructor.
    MetisDivideHeterogeneousInputProcess(MetisDivideHeterogeneousInputProcess const& rOther):
        MetisDivideInputToPartitionsProcess(rOther.mrIO,rOther.mNumberOfPartitions,rOther.mDimension),
        mSynchronizeConditions(rOther.mSynchronizeConditions),
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
    void Execute() override
    {
        // Read nodal graph from input

        IO::ConnectivitiesContainerType KratosFormatNodeConnectivities;

        SizeType NumNodes = BaseType::mrIO.ReadNodalGraph(KratosFormatNodeConnectivities);

        SizeType NumNodesInMesh = BaseType::mrIO.ReadNodesNumber();
        if (NumNodes != NumNodesInMesh)
            KRATOS_ERROR << "Invalid mesh: number of connected nodes = " << NumNodes
                         << ", number of mesh nodes = " << NumNodesInMesh << "."
                         << std::endl;

        // Write connectivity data in CSR format
        idxtype* NodeIndices = 0;
        idxtype* NodeConnectivities = 0;

        ConvertKratosToCSRFormat(KratosFormatNodeConnectivities, &NodeIndices, &NodeConnectivities);

        std::vector<idxtype> NodePartition;
        PartitionNodes(NumNodes,NodeIndices,NodeConnectivities,NodePartition);

        // Free some memory we no longer need
        delete [] NodeIndices;
        delete [] NodeConnectivities;

        // Partition elements
        IO::ConnectivitiesContainerType ElementConnectivities;
        SizeType NumElements =  BaseType::mrIO.ReadElementsConnectivities(ElementConnectivities);
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
        SizeType NumConditions = BaseType::mrIO.ReadConditionsConnectivities(ConditionConnectivities);
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
        CalculateDomainsGraph(DomainGraph,NumElements,ElementConnectivities,NodePartition,ElementPartition);
        CalculateDomainsGraph(DomainGraph,NumConditions,ConditionConnectivities,NodePartition,ConditionPartition);

        int NumColors;
        GraphType ColoredDomainGraph;
        GraphColoringProcess(mNumberOfPartitions,DomainGraph,ColoredDomainGraph,NumColors).Execute();

        if (mVerbosity > 0)
        {
            KRATOS_WATCH(NumColors);
        }

	if (mVerbosity > 2)
	{
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
    std::string Info() const override
    {
        return "MetisDivideHeterogeneousInputProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MetisDivideHeterogeneousInputProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
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

protected:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    bool mSynchronizeConditions;

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
                       idxtype* NodeIndices,
                       idxtype* NodeConnectivities,
                       std::vector<idxtype>& rNodePartition)
    {
        idxtype n = static_cast<idxtype>(NumNodes);

        idxtype nparts = static_cast<idxtype>(BaseType::mNumberOfPartitions);
        idxtype edgecut;
        rNodePartition.resize(NumNodes);


#ifndef KRATOS_USE_METIS_5
        idxtype wgtflag = 0; // Graph is not weighted
        idxtype numflag = 0; // Nodes are numbered from 0 (C style)
        idxtype options[5]; // options array
        options[0] = 0; // use default options

        //old version
        METIS_PartGraphKway(&n,NodeIndices,NodeConnectivities,NULL,NULL,&wgtflag,&numflag,&nparts,&options[0],&edgecut,&rNodePartition[0]);
#else
        idx_t ncon = 1; //The number of balancing constraints. It should be at least 1.

        idx_t options[METIS_NOPTIONS];
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
#endif



        // Debug: print partition
        PrintDebugData("Node Partition",rNodePartition);

        return edgecut;
    }

    /// Use the nodal partition data to assign elements or conditions to a partition.
    void PartitionMesh(std::vector<idxtype> const& NodePartition,
                       const IO::ConnectivitiesContainerType& rElemConnectivities,
                       std::vector<idxtype>& rElemPartition)
    {
        SizeType NumElements = rElemConnectivities.size();
        std::vector<int> PartitionWeights(BaseType::mNumberOfPartitions,0);

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
        int MaxWeight = 1.03 * NumElements / BaseType::mNumberOfPartitions;
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

    /// Partition the elements such that boundary elements are always assigned the majority partition.
    void PartitionElementsSynchronous(std::vector<idxtype> const& NodePartition,
                       const IO::ConnectivitiesContainerType& rElemConnectivities,
                       std::vector<idxtype>& rElemPartition)
    {
        SizeType NumElements = rElemConnectivities.size();
        std::vector<int> PartitionWeights(BaseType::mNumberOfPartitions,0);

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
        //int MaxWeight = 1.03 * NumElements / BaseType::mNumberOfPartitions;
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
                    PartitionWeights[MajorityPartition]++;
                }
            }

            // Advance to next element in connectivities array
            itElem++;
        }

        PrintDebugData("Element Partition",rElemPartition);

    }

    /// Partition the conditions such that the condition is assigned the same partition as its parent element.
    void PartitionConditionsSynchronous(const std::vector<idxtype>& rNodePartition,
			     const std::vector<idxtype>& rElemPartition,
			     const IO::ConnectivitiesContainerType& rCondConnectivities,
			     const IO::ConnectivitiesContainerType& rElemConnectivities,
			     std::vector<idxtype>& rCondPartition)
    {
      SizeType NumElements = rElemConnectivities.size();
      SizeType NumConditions = rCondConnectivities.size();
      std::vector<int> PartitionWeights(BaseType::mNumberOfPartitions,0);

      // initialize CondPartition
      rCondPartition.resize(NumConditions,-1);

      // make sorted element connectivities array
      IO::ConnectivitiesContainerType ElementsSorted(rElemConnectivities);
      for (SizeType i=0; i<NumElements; i++)
          std::sort(ElementsSorted[i].begin(), ElementsSorted[i].end());

      // Conditions where all nodes belong to the same partition always go to that partition
      IO::ConnectivitiesContainerType::const_iterator itCond = rCondConnectivities.begin();
      for (std::vector<idxtype>::iterator itPart = rCondPartition.begin(); itPart != rCondPartition.end(); itPart++)
      {
          int MyPartition = rNodePartition[ (*itCond)[0] - 1 ]; // Node Ids start from 1
          SizeType NeighbourNodes = 1; // Nodes in the same partition
          for (std::vector<SizeType>::const_iterator itNode = itCond->begin()+1; itNode != itCond->end(); ++itNode)
          {
              if ( rNodePartition[ *itNode - 1 ] == MyPartition )
                  ++NeighbourNodes;
              else
                  break;
          }

          if ( NeighbourNodes == itCond->size() )
          {
              *itPart = MyPartition;
              PartitionWeights[MyPartition]++;
          }

          // Advance to next condition in connectivities array
          itCond++;
      }
      // Now distribute boundary conditions
      itCond = rCondConnectivities.begin();
      //int MaxWeight = 1.03 * NumConditions / BaseType::mNumberOfPartitions;
      for (std::vector<idxtype>::iterator itPart = rCondPartition.begin(); itPart != rCondPartition.end(); itPart++)
      {
          if (*itPart == -1) // If condition is still unassigned
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
                  PartitionWeights[MajorityPartition]++;
              }

              // ensure conditions sharing nodes with an element have same partition as the element
              IO::ConnectivitiesContainerType::value_type tmp(*itCond);
              std::sort(tmp.begin(), tmp.end());

              for (SizeType i=0; i<NumElements; i++)
              {
                  if ( std::includes(ElementsSorted[i].begin(), ElementsSorted[i].end(), tmp.begin(), tmp.end()) )
                  {
                      *itPart = rElemPartition[i];
                      break;
                  }
              }
          }

          // Advance to next condition in connectivities array
          itCond++;
      }

      PrintDebugData("Condition Partition",rCondPartition);

    }

    void RedistributeHangingNodes(
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
                        const std::vector<idxtype>& rPartitionData)
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
                std::cout << count << " objects." << std::endl;
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
