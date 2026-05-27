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

// System includes
#include <fstream>
#include <sstream>

// Project includes
#include "custom_processes/kahip_divide_heterogeneous_input_in_memory_process.h"

namespace Kratos {

// ─── Constructors ─────────────────────────────────────────────────────────────

KaHIPDivideHeterogeneousInputInMemoryProcess::KaHIPDivideHeterogeneousInputInMemoryProcess(
    IO& rIO,
    ModelPartIO& rSerialIO,
    const DataCommunicator& rDataComm,
    Parameters rSettings,
    bool SynchronizeConditions)
    : BaseType(rIO, static_cast<SizeType>(rDataComm.Size()), rSettings, SynchronizeConditions),
      mrSerialIO(rSerialIO),
      mrDataComm(rDataComm)
{
    KRATOS_ERROR_IF_NOT(mrDataComm.IsDistributed())
        << "KaHIPDivideHeterogeneousInputInMemoryProcess: "
        << "DataCommunicator must be distributed (MPI). "
        << "For single-process runs use KaHIPDivideHeterogeneousInputProcess." << std::endl;
}

KaHIPDivideHeterogeneousInputInMemoryProcess::KaHIPDivideHeterogeneousInputInMemoryProcess(
    IO& rIO,
    ModelPartIO& rSerialIO,
    const DataCommunicator& rDataComm,
    int Dimension,
    int Verbosity,
    bool SynchronizeConditions)
    : BaseType(rIO, static_cast<SizeType>(rDataComm.Size()), Dimension, Verbosity, SynchronizeConditions),
      mrSerialIO(rSerialIO),
      mrDataComm(rDataComm)
{
    KRATOS_ERROR_IF_NOT(mrDataComm.IsDistributed())
        << "KaHIPDivideHeterogeneousInputInMemoryProcess: "
        << "DataCommunicator must be distributed (MPI). "
        << "For single-process runs use KaHIPDivideHeterogeneousInputProcess." << std::endl;
}

// ─── Execute ──────────────────────────────────────────────────────────────────

void KaHIPDivideHeterogeneousInputInMemoryProcess::Execute()
{
    KRATOS_TRY

    const int mpi_rank = mrDataComm.Rank();
    const int mpi_size = mrDataComm.Size();

    std::vector<std::vector<char>> send_buffer;

    // ── Rank 0: partition and serialise each per-rank model part ────────────
    if (mpi_rank == 0) {
        PartitioningInfo part_info;
        ExecutePartitioning(part_info);

        // Write each partition into its own stringstream
        std::vector<Kratos::shared_ptr<std::iostream>> streams(mpi_size);
        std::vector<std::stringbuf> stringbufs(mpi_size);
        for (int i = 0; i < mpi_size; ++i) {
            streams[i] = Kratos::make_shared<std::iostream>(&stringbufs[i]);
        }

        mrIO.DivideInputToPartitions(
            streams.data(),
            mNumberOfPartitions,
            part_info);

        // Convert stringbufs to char vectors for Scatterv
        send_buffer.resize(mpi_size);
        for (int i = 0; i < mpi_size; ++i) {
            const std::string str = stringbufs[i].str();
            send_buffer[i].reserve(str.size());
            std::copy(str.begin(), str.end(), std::back_inserter(send_buffer[i]));
        }
    }

    // ── All ranks: receive the local partition buffer ─────────────────────
    const auto recv_buffer = mrDataComm.Scatterv(send_buffer, 0);

    // Wrap received data in a stream and hand it to the serial IO
    auto p_local_stream = Kratos::make_shared<std::iostream>(new std::stringbuf());
    p_local_stream->write(recv_buffer.data(), recv_buffer.size());

    if (mVerbosity > 1) {
        std::ofstream debug_file(
            "KaHIPDivideHeterogeneousInputInMemoryProcess_debug_modelpart_"
            + std::to_string(mpi_rank) + ".mdpa");
        debug_file.write(recv_buffer.data(), static_cast<std::streamsize>(recv_buffer.size()));
    }

    // Redirect the serial IO to read from the in-memory buffer
    mrSerialIO.SwapStreamSource(p_local_stream);

    KRATOS_CATCH("")
}

} // namespace Kratos
