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
    MPI_Comm this_comm = MPIDataCommunicator::GetMPICommunicator(mrDataComm);

    int mpi_rank = mrDataComm.Rank();
    int mpi_size = mrDataComm.Size();

    int * msgSendSize = new int[mpi_size];
    int * msgRecvSize = new int[mpi_size];

    const char ** mpi_send_buffer = new const char * [mpi_size];
    char ** mpi_recv_buffer = new char * [mpi_size];
    std::string * str = new std::string[mpi_size];

    // Set size
    for(int i = 0; i < mpi_size; i++) {
        msgSendSize[i] = 0;
        msgRecvSize[i] = 0;
    }

    // Transfer Streams
    Kratos::shared_ptr<std::iostream> * streams = new Kratos::shared_ptr<std::iostream>[mpi_size];
    std::stringbuf * stringbufs = new std::stringbuf[mpi_size];

    for(auto i = 0; i < mpi_size; i++) {
        streams[i] = Kratos::shared_ptr<std::iostream>(new std::iostream(&stringbufs[i]));
    }

    // Main process calculates the partitions and writes the result into temporal streams
    if(mpi_rank == 0) {

        auto part_info(Kratos::make_shared<PartitioningInfo>());

        ExecutePartitioning(*part_info);

        // Write partitions to streams
        mrIO.DivideInputToPartitions(
            streams,
            mNumberOfPartitions,
            part_info->Graph,
            part_info->NodesPartitions,
            part_info->ElementsPartitions,
            part_info->ConditionsPartitions,
            part_info->NodesAllPartitions,
            part_info->ElementsAllPartitions,
            part_info->ConditionsAllPartitions);
    }

    // Calculate the message and prepare the buffers
    if(mpi_rank == 0) {
        for(auto i = 0; i < mpi_size; i++) {
            str[i] = stringbufs[i].str();
            msgSendSize[i] = str[i].size();
            mpi_send_buffer[i] = str[i].c_str();
        }
    }

    // Send the message size to all processes
    MPI_Scatter(msgSendSize,1,MPI_INT,&msgRecvSize[mpi_rank],1,MPI_INT,0,this_comm);

    // Calculate the number of events:
    auto NumberOfCommunicationEvents = 1 + mpi_size * !mpi_rank;
    auto NumberOfCommunicationEventsIndex = 0;

    // Prepare the communication events
    MPI_Request * reqs = new MPI_Request[NumberOfCommunicationEvents];
    MPI_Status * stats = new MPI_Status[NumberOfCommunicationEvents];

    // Set up all receive and send events
    if( mpi_rank == 0) {
        for(auto i = 0; i < mpi_size; i++) {
            char* aux_char = const_cast<char*>(mpi_send_buffer[i]);
            MPI_Isend(aux_char,msgSendSize[i],MPI_CHAR,i,0,this_comm,&reqs[NumberOfCommunicationEventsIndex++]);
        }
    }

    // Recieve the buffers
    mpi_recv_buffer[mpi_rank] = (char *)malloc(sizeof(char) * msgRecvSize[mpi_rank]);
    MPI_Irecv(mpi_recv_buffer[mpi_rank],msgRecvSize[mpi_rank],MPI_CHAR,0,0,this_comm,&reqs[NumberOfCommunicationEventsIndex++]);

    // Wait untill all communications finish
    if( MPI_Waitall(NumberOfCommunicationEvents, reqs, stats) != MPI_SUCCESS ) {
        KRATOS_ERROR << "Error in metis_partition_mem" << std::endl;
    }

    mrDataComm.Barrier();

    if(mpi_rank != 0) {
        streams[mpi_rank]->write(mpi_recv_buffer[mpi_rank], msgRecvSize[mpi_rank]);
    }

    if (mVerbosity > 1) {
        std::ofstream debug_ofstream("MetisDivideHeterogeneousInputInMemoryProcess_debug_modelpart_"+std::to_string(mpi_rank)+".mdpa");
        debug_ofstream << stringbufs[mpi_rank].str() << std::endl;
    }

    // TODO: Try to come up with a better way to change the buffer.
    mrSerialIO.SwapStreamSource(streams[mpi_rank]);

    // Free buffers
    free(mpi_recv_buffer[mpi_rank]);

    delete [] reqs;
    delete [] stats;

    delete [] mpi_recv_buffer;
    delete [] mpi_send_buffer;
    delete [] str;

    delete [] msgSendSize;
    delete [] msgRecvSize;
}

} // namespace Kratos
