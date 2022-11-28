#pragma once

// External includes

// Project includes
#include "metis_divide_input_to_partitions_process.h"

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
class KRATOS_API(METIS_APPLICATION) MetisDivideHeterogeneousInputProcess : public MetisDivideInputToPartitionsProcess
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MetisDivideHeterogeneousInputProcess
    KRATOS_CLASS_POINTER_DEFINITION(MetisDivideHeterogeneousInputProcess);

    using SizeType = IO::SizeType;
    using GraphType = IO::GraphType;
    using PartitionIndicesType = IO::PartitionIndicesType;
    using PartitionIndicesContainerType = IO::PartitionIndicesContainerType;
    using idxtype = idx_t; // from metis

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
    void Execute() override;

    virtual void GetNodesPartitions(std::vector<idxtype> &rNodePartition, SizeType &rNumNodes);

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
    ///@name Protected LifeCycle
    ///@{

    struct PartitioningInfo
    {
        GraphType Graph;
        PartitionIndicesType NodesPartitions; // partition where the Node is local
        PartitionIndicesType ElementsPartitions; // partition where the Element is local
        PartitionIndicesType ConditionsPartitions; // partition where the Condition is local
        PartitionIndicesContainerType NodesAllPartitions; // partitions, in which the Node is present (local & ghost)
        PartitionIndicesContainerType ElementsAllPartitions; // partitions, in which the Element is present (local & ghost)
        PartitionIndicesContainerType ConditionsAllPartitions; // partitions, in which the Condition is present (local & ghost)
    };

    ///@}
    ///@name Member Variables
    ///@{

    bool mSynchronizeConditions;

    int mVerbosity;
    int mNumNodes;

    std::vector<std::unordered_set<std::size_t>> mNodeConnectivities;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void ExecutePartitioning(PartitioningInfo& rPartitioningInfo);

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
                       std::vector<idxtype>& rNodePartition);

    /// Use the nodal partition data to assign elements or conditions to a partition.
    void PartitionMesh(std::vector<idxtype> const& NodePartition,
                       const IO::ConnectivitiesContainerType& rElemConnectivities,
                       std::vector<idxtype>& rElemPartition);

    /// Partition the elements such that boundary elements are always assigned the majority partition.
    void PartitionElementsSynchronous(std::vector<idxtype> const& NodePartition,
                       const IO::ConnectivitiesContainerType& rElemConnectivities,
                       std::vector<idxtype>& rElemPartition);

    /// Partition the conditions such that the condition is assigned the same partition as its parent element.
    void PartitionConditionsSynchronous(const std::vector<idxtype>& rNodePartition,
			     const std::vector<idxtype>& rElemPartition,
			     const IO::ConnectivitiesContainerType& rCondConnectivities,
			     const IO::ConnectivitiesContainerType& rElemConnectivities,
			     std::vector<idxtype>& rCondPartition);

    void RedistributeHangingNodes(
            std::vector<idxtype>& rNodePartition,
            std::vector<idxtype> const& rElementPartition,
            const IO::ConnectivitiesContainerType& rElementConnectivities,
            std::vector<idxtype> const& rConditionPartition,
            const IO::ConnectivitiesContainerType& rConditionConnectivities);

    SizeType FindMax(SizeType NumTerms, const std::vector<int>& rVect);

    void PrintDebugData(const std::string& rLabel,
                        const std::vector<idxtype>& rPartitionData);

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
