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

#pragma once

// Project includes
#include "includes/data_communicator.h"
#include "includes/model_part_io.h"
#include "custom_processes/kahip_divide_heterogeneous_input_process.h"

namespace Kratos
{

///@addtogroup KaHIPApplication
///@{

///@name Kratos Classes
///@{

/**
 * @class KaHIPDivideHeterogeneousInputInMemoryProcess
 * @ingroup KaHIPApplication
 * @brief MPI-parallel in-memory mesh partitioning using KaHIP.
 * @details Drop-in replacement for @c MetisDivideHeterogeneousInputInMemoryProcess.
 *
 *          Partitioning is performed exclusively on MPI rank 0 using the serial
 *          @c KaHIPDivideHeterogeneousInputProcess base. After partitioning, the
 *          per-rank .mdpa data is serialised into strings and distributed to all
 *          ranks via @c DataCommunicator::Scatterv(). Each rank then redirects its
 *          @c ModelPartIO to read from the received in-memory buffer by calling
 *          @c ModelPartIO::SwapStreamSource().
 *
 *          This process is used when @c partition_in_memory = @c true, avoiding
 *          the I/O overhead of writing intermediate per-rank files to disk.
 *
 *          **Compatibility**: constructor signatures mirror
 *          @c MetisDivideHeterogeneousInputInMemoryProcess so that existing Python
 *          scripts only need to change the import and class name.
 *
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(KAHIP_APPLICATION) KaHIPDivideHeterogeneousInputInMemoryProcess
    : public KaHIPDivideHeterogeneousInputProcess
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(KaHIPDivideHeterogeneousInputInMemoryProcess);

    using BaseType = KaHIPDivideHeterogeneousInputProcess;
    using BaseType::SizeType;
    using BaseType::idxtype;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Parameters-based constructor.
     * @param rIO                IO object to read mesh connectivity from (used on rank 0).
     * @param rSerialIO          ModelPartIO that will be redirected after scatter to hold
     *                           the local partition's .mdpa data.
     * @param rDataComm          Distributed MPI data communicator.
     * @param rSettings          KaHIP partitioner settings (see KaHIPPartitioner::GetDefaultParameters).
     * @param SynchronizeConditions  Co-locate conditions with their parent element partition.
     */
    KaHIPDivideHeterogeneousInputInMemoryProcess(
        IO& rIO,
        ModelPartIO& rSerialIO,
        const DataCommunicator& rDataComm,
        Parameters rSettings,
        bool SynchronizeConditions = false);

    /**
     * @brief Compatibility constructor matching MetisDivideHeterogeneousInputInMemoryProcess.
     * @details Uses ECO mode with default imbalance 0.03 and a single trial.
     */
    KaHIPDivideHeterogeneousInputInMemoryProcess(
        IO& rIO,
        ModelPartIO& rSerialIO,
        const DataCommunicator& rDataComm,
        int Dimension = 3,
        int Verbosity = 0,
        bool SynchronizeConditions = false);

    /// Destructor.
    ~KaHIPDivideHeterogeneousInputInMemoryProcess() override = default;

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
     * @brief Partition on rank 0, then Scatterv per-rank streams to all ranks.
     * @details
     *  - Rank 0: calls @c ExecutePartitioning(), serialises each partition into
     *    a @c stringstream, builds a send-buffer for @c Scatterv.
     *  - All ranks: receive their buffer and call
     *    @c mrSerialIO.SwapStreamSource() so that subsequent @c ReadModelPart
     *    calls on @c mrSerialIO read from the in-memory buffer rather than disk.
     */
    void Execute() override;

    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const override
    {
        return "KaHIPDivideHeterogeneousInputInMemoryProcess";
    }

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    void PrintData(std::ostream& rOStream) const override {}

    ///@}

private:
    ///@name Member Variables
    ///@{

    ModelPartIO& mrSerialIO;
    const DataCommunicator& mrDataComm;

    ///@}

}; // class KaHIPDivideHeterogeneousInputInMemoryProcess

///@}

} // namespace Kratos
