//     __ __      __  __________  ___                ___            __  _           
//    / //_/___ _/ / / /  _/ __ \/   |  ____  ____  / (_)________ _/ /_(_)___  ____ 
//   / ,< / __ `/ /_/ // // /_/ / /| | / __ \/ __ \/ / / ___/ __ `/ __/ / __ \/ __ \
//  / /| / /_/ / __  // // ____/ ___ |/ /_/ / /_/ / / / /__/ /_/ / /_/ / / /_/ / / /
// /_/ |_\__,_/_/ /_/___/_/   /_/  |_/ .___/ .___/_/_/\___/\__,_/\__/_/\____/_/ /_/ 
//                                  /_/   /_/                                       
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes
#include "kaHIP_interface.h"

// Project includes
#include "includes/io.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "custom_utilities/kahip_partitioner.h"

namespace Kratos
{

///@addtogroup KaHIPApplication
///@{

///@name Kratos Classes
///@{

/**
 * @class KaHIPDivideHeterogeneousInputProcess
 * @ingroup KaHIPApplication
 * @brief Partitions a heterogeneous Kratos mesh read through an @c IO interface using KaHIP.
 * @details This process is the KaHIP drop-in replacement for
 *          @c MetisDivideHeterogeneousInputProcess. It:
 *
 *          1. Reads the nodal connectivity graph from the provided @c IO object via
 *             @c IO::ReadNodalGraph().
 *          2. Converts the graph to CSR format using KaHIPCSRConverter.
 *          3. Calls @c kaffpa() via KaHIPPartitioner to partition the nodal graph into
 *             @p NumberOfPartitions blocks.
 *          4. Assigns elements, conditions, geometries, and master-slave constraints to
 *             partitions based on the majority-partition of their nodes.
 *          5. Optionally synchronises conditions so that boundary conditions are always
 *             co-located with their parent element (when @p SynchronizeConditions is
 *             @c true).
 *          6. Writes the partitioned model part files by calling
 *             @c IO::DivideInputToPartitions().
 *
 *          **Compatibility note**: the constructor signatures mirror
 *          @c MetisDivideHeterogeneousInputProcess so that existing Python scripts
 *          only need to change the import and class name.
 *
 *          **KaHIP preconfiguration** can be set via the @p rSettings @c Parameters
 *          overload. When using the legacy integer-based constructors, the ECO mode is
 *          used by default.
 *
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(KAHIP_APPLICATION) KaHIPDivideHeterogeneousInputProcess 
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KaHIPDivideHeterogeneousInputProcess
    KRATOS_CLASS_POINTER_DEFINITION(KaHIPDivideHeterogeneousInputProcess);

    /// Some IO definitions
    using SizeType                       = IO::SizeType;
    using GraphType                      = IO::GraphType;
    using PartitioningInfo               = IO::PartitioningInfo;
    using PartitionIndicesType           = IO::PartitionIndicesType;
    using PartitionIndicesContainerType  = IO::PartitionIndicesContainerType;

    /// KaHIP index type (int32_t by default; int64_t with 64-bit build)
    using idxtype = kahip_idx;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Main constructor accepting KaHIP configuration via Parameters.
     * @param rIO                  IO object to read mesh connectivity from and write partitions to.
     * @param NumberOfPartitions   Number of blocks k to partition the mesh into.
     * @param rSettings            KaHIP partitioner settings (see KaHIPPartitioner::GetDefaultParameters).
     *                             Additional keys accepted here:
     *                             - "dimension" (int, default 3): spatial dimension (used by IO)
     *                             - "verbosity" (int, default 0): 0 = silent, 1 = info, 2 = debug
     *                             - "synchronize_conditions" (bool, default false): co-locate
     *                               conditions with their parent element
     * @param SynchronizeConditions  Legacy shortcut to enable condition synchronisation.
     */
    KaHIPDivideHeterogeneousInputProcess(
        IO& rIO,
        SizeType NumberOfPartitions,
        Parameters rSettings,
        bool SynchronizeConditions = false);

    /**
     * @brief Compatibility constructor matching MetisDivideHeterogeneousInputProcess signature.
     * @details Uses ECO mode with default imbalance 0.03 and a single trial.
     */
    KaHIPDivideHeterogeneousInputProcess(
        IO& rIO,
        SizeType NumberOfPartitions,
        int Dimension = 3,
        int Verbosity = 0,
        bool SynchronizeConditions = false);

    /// Destructor.
    ~KaHIPDivideHeterogeneousInputProcess() override = default;

    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Partition the mesh and write per-rank .mdpa files.
     * @details Performs the complete partitioning workflow described in the class documentation.
     */
    void Execute() override;

    /**
     * @brief Partition the nodal graph and store results in @p rNodePartition.
     * @details Reads the nodal graph from @c mrIO, converts to CSR, calls @c kaffpa(),
     *          and fills @p rNodePartition with the block index for each node (0-indexed).
     * @param rNodePartition  [out] Block assignment for each node (size = num_nodes)
     * @param rNumNodes       [out] Total number of nodes in the mesh
     */
    virtual void GetNodesPartitions(
        std::vector<idxtype>& rNodePartition,
        SizeType& rNumNodes
        );

    ///@}
    ///@name Input and output
    ///@{

    /**
     * @brief Retrieves the name of the class
     * @return Name of the class string 
     */
    std::string Info() const override
    {
        return "KaHIPDivideHeterogeneousInputProcess";
    }

    /**
     * @brief Prints the info to stream
     * @param rOStream The stream considered
     */
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /**
     * @brief Prints the data to stream
     * @param rOStream The stream considered
     */
    void PrintData(std::ostream& rOStream) const override 
    {
    }

    ///@}

protected:
    ///@name Member Variables
    ///@{

    IO&          mrIO;                     /// The IO instance considered
    SizeType     mNumberOfPartitions;      /// The number of partitions
    bool         mSynchronizeConditions;   /// The bool to know if the conditiosn are synced
    int          mVerbosity;               /// Verbosity level
    int          mNumNodes;                /// The number of nodes

    /// Adjacency sets built while reading the nodal graph (cached for subclass access)
    std::vector<std::unordered_set<std::size_t>> mNodeConnectivities;

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief Execute the partitioning logic and fill @p rPartitioningInfo.
     * @details This method is separated from Execute() to allow subclasses (e.g. the
     *          in-memory MPI variant) to reuse it without writing files.
     */
    void ExecutePartitioning(PartitioningInfo& rPartitioningInfo);

    ///@}

private:
    ///@name Member Variables
    ///@{

    KaHIPPartitioner mPartitioner; /// The KaHIP partitioner used

    ///@}
    ///@name Private Operations
    ///@{

    /// Assign elements/conditions/geometries to the partition that contains the majority
    /// of their nodes (ties broken by lowest partition index).
    void PartitionMesh(
        const std::vector<idxtype>& rNodePartition,
        const IO::ConnectivitiesContainerType& rConnectivities,
        std::vector<idxtype>& rEntityPartition);

    /// Variant that ensures geometries share the partition of their majority node.
    void PartitionGeometriesSynchronous(
        const std::vector<idxtype>& rNodePartition,
        const IO::ConnectivitiesContainerType& rGeometryConnectivities,
        std::vector<idxtype>& rGeometryPartition);

    /// Variant that ensures elements sharing nodes on partition boundaries are kept
    /// in the partition that contains the most nodes.
    void PartitionElementsSynchronous(
        const std::vector<idxtype>& rNodePartition,
        const IO::ConnectivitiesContainerType& rElemConnectivities,
        std::vector<idxtype>& rElemPartition);

    /// Assign each condition to the same partition as its parent element.
    void PartitionConditionsSynchronous(
        const std::vector<idxtype>& rNodePartition,
        const std::vector<idxtype>& rElemPartition,
        const IO::ConnectivitiesContainerType& rCondConnectivities,
        const IO::ConnectivitiesContainerType& rElemConnectivities,
        std::vector<idxtype>& rCondPartition);

    /// Assign constraints to the partition of their majority node.
    void PartitionConstraintsSynchronous(
        const std::vector<idxtype>& rNodePartition,
        const IO::ConnectivitiesContainerType& rConstraintConnectivities,
        std::vector<idxtype>& rConstraintPartition);

    /// Move hanging nodes (nodes with no local element in their partition) to a
    /// partition that does have an element using them.
    void RedistributeHangingNodes(
        std::vector<idxtype>& rNodePartition,
        const std::vector<idxtype>& rGeometryPartition,
        const IO::ConnectivitiesContainerType& rGeometryConnectivities,
        const std::vector<idxtype>& rElementPartition,
        const IO::ConnectivitiesContainerType& rElementConnectivities,
        const std::vector<idxtype>& rConditionPartition,
        const IO::ConnectivitiesContainerType& rConditionConnectivities,
        const std::vector<idxtype>& rConstraintPartition,
        const IO::ConnectivitiesContainerType& rConstraintConnectivities);

    SizeType FindMax(SizeType NumTerms, const std::vector<int>& rVect);

    void PrintDebugData(
        const std::string& rLabel,
        const std::vector<idxtype>& rPartitionData);

    ///@}

}; // class KaHIPDivideHeterogeneousInputProcess

///@}

} // namespace Kratos
