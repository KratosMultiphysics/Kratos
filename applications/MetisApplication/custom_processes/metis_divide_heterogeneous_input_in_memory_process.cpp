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
#include "mpi/includes/mpi_data_communicator.h"
#include "metis_divide_heterogeneous_input_in_memory_process.h"


namespace Kratos {

void MetisDivideHeterogeneousInputInMemoryProcess::Execute()
{
    KRATOS_TRY

    MPI_Comm this_comm = MPIDataCommunicator::GetMPICommunicator(mrDataComm);

    int mpi_rank = mrDataComm.Rank();
    int mpi_size = mrDataComm.Size();

    // const char ** mpi_send_buffer = new const char * [mpi_size];
    // std::string * str = new std::string[mpi_size];



    using BufferType = std::vector<char>;

    std::vector<BufferType> send_buffer;

    // Main process calculates the partitions and writes the result into temporal streams
    if(mpi_rank == 0) {

        auto part_info(Kratos::make_shared<PartitioningInfo>());

        ExecutePartitioning(*part_info);

        // Transfer Streams
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

    KRATOS_WATCH("1111");

    BufferType recv_buffer = mrDataComm.Scatterv(send_buffer, 0);
    KRATOS_WATCH("222");

    // // Calculate the message and prepare the buffers
    // std::vector<int> send_sizes(mpi_size);
    // if(mpi_rank == 0) {
    //     for(auto i = 0; i < mpi_size; i++) {
    //         str[i] = stringbufs[i].str();
    //         send_sizes[i] = str[i].size();
    //         mpi_send_buffer[i] = str[i].c_str();
    //     }
    // }

    // // Send the message size to all processes
    // std::vector<int> vec_recv_size(1);
    // mrDataComm.Scatter(send_sizes, vec_recv_size, 0);
    // const int recv_size = vec_recv_size[0];

    // // Calculate the number of events:
    // auto NumberOfCommunicationEvents = 1 + mpi_size * !mpi_rank;
    // auto NumberOfCommunicationEventsIndex = 0;

    // // Prepare the communication events
    // MPI_Request * reqs = new MPI_Request[NumberOfCommunicationEvents];
    // MPI_Status * stats = new MPI_Status[NumberOfCommunicationEvents];

    // // Set up all receive and send events
    // if( mpi_rank == 0) {
    //     for(auto i = 0; i < mpi_size; i++) {
    //         char* aux_char = const_cast<char*>(mpi_send_buffer[i]);
    //         MPI_Isend(aux_char,send_sizes[i],MPI_CHAR,i,0,this_comm,&reqs[NumberOfCommunicationEventsIndex++]);
    //     }
    // }

    // // Receive the buffers
    // std::vector<char> recv_buffer(recv_size);
    // MPI_Irecv(recv_buffer.data(), recv_size,MPI_CHAR,0,0,this_comm,&reqs[NumberOfCommunicationEventsIndex++]);

    // // Wait untill all communications finish
    // if( MPI_Waitall(NumberOfCommunicationEvents, reqs, stats) != MPI_SUCCESS ) {
    //     KRATOS_ERROR << "Error in metis_partition_mem" << std::endl;
    // }

    // mrDataComm.Barrier();

    // if(mpi_rank != 0) {
        // streams[mpi_rank]->write(recv_buffer.data(), recv_buffer.size());
    // }
    KRATOS_WATCH("333");

    // if (mVerbosity > 1) {
    //     std::ofstream debug_ofstream("MetisDivideHeterogeneousInputInMemoryProcess_debug_modelpart_"+std::to_string(mpi_rank)+".mdpa");
    //     debug_ofstream << stringbufs[mpi_rank].str() << std::endl;
    // }
    KRATOS_WATCH("444");

    auto p_local_stream(Kratos::make_shared<std::iostream>(new std::stringbuf()));
    p_local_stream->write(recv_buffer.data(), recv_buffer.size());

    // TODO: Try to come up with a better way to change the buffer.
    mrSerialIO.SwapStreamSource(p_local_stream);
    KRATOS_WATCH("555");

    // Free buffers
    // delete [] reqs;
    // delete [] stats;

    // delete [] mpi_send_buffer;
    // delete [] str;
    KRATOS_WATCH("666");

    KRATOS_CATCH("")
}

} // namespace Kratos
