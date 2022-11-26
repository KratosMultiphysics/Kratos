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


// Project includes
#include "metis_divide_heterogeneous_input_in_memory_process.h"


namespace Kratos {

void MetisDivideHeterogeneousInputInMemoryProcess::Execute()
{
    KRATOS_TRY

    const int mpi_rank = mrDataComm.Rank();
    const int mpi_size = mrDataComm.Size();

    std::vector<std::vector<char>> send_buffer;

    // Main process calculates the partitions and writes the result into temporal streams
    if(mpi_rank == 0) {

        auto part_info(Kratos::make_shared<PartitioningInfo>());

        ExecutePartitioning(*part_info);

        // prepare partitioning streams
        std::vector<Kratos::shared_ptr<std::iostream>> streams(mpi_size);
        std::vector<std::stringbuf> stringbufs(mpi_size);

        for(int i = 0; i < mpi_size; i++) {
            streams[i] = Kratos::make_shared<std::iostream>(&stringbufs[i]);
        }

        // Write partitions to streams
        mrIO.DivideInputToPartitions(
            streams.data(),
            mNumberOfPartitions,
            part_info->Graph,
            part_info->NodesPartitions,
            part_info->ElementsPartitions,
            part_info->ConditionsPartitions,
            part_info->NodesAllPartitions,
            part_info->ElementsAllPartitions,
            part_info->ConditionsAllPartitions);

        // prepare buffer to scatter to other partitions
        send_buffer.resize(mpi_size);
        for(int i=0; i<mpi_size; ++i) {
            const std::string str = stringbufs[i].str();
            send_buffer[i].reserve(str.size());
            std::copy(str.begin(), str.end(), std::back_inserter(send_buffer[i]));
        }
    }

    const auto recv_buffer = mrDataComm.Scatterv(send_buffer, 0);

    auto p_local_stream(Kratos::make_shared<std::iostream>(new std::stringbuf()));
    p_local_stream->write(recv_buffer.data(), recv_buffer.size());

    if (mVerbosity > 1) {
        std::ofstream debug_ofstream("MetisDivideHeterogeneousInputInMemoryProcess_debug_modelpart_"+std::to_string(mpi_rank)+".mdpa");
        debug_ofstream << recv_buffer << std::endl;
    }

    // TODO: Try to come up with a better way to change the buffer.
    mrSerialIO.SwapStreamSource(p_local_stream);

    KRATOS_CATCH("")
}

} // namespace Kratos
