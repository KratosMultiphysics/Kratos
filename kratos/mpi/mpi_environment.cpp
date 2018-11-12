//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Jordi Cotela
//

// System includes
#include <mpi.h>

// Project includes
#include "input_output/logger.h"
#include "includes/parallel_environment.h"

// Module includes
#include "mpi_environment.h"
#include "mpi/includes/mpi_data_communicator.h"

namespace Kratos {

void MPIEnvironment::Initialize()
{
    if (!MPIEnvironment::IsInitialized())
    {
        // Note from openmpi docs: https://www.open-mpi.org/doc/v2.0/man3/MPI_Init.3.php
        // "Open MPI accepts the C/C++ argc and argv arguments to main, but neither modifies, interprets, nor distributes them"
        // I am passing empty arguments for simplicity.
        int argc = 0;
        char** argv = nullptr;

        #if MPI_VERSION < 2
        MPI_Init(&argc, &argv);

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) KRATOS_DETAIL("MPIEnvironment") << "MPI version older than 2, no support for MPI_THREAD_MULTIPLE is provided." << std::endl;

        #else
        int provided;
        MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

        if (provided < MPI_THREAD_MULTIPLE)
        {
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            if (rank == 0) KRATOS_DETAIL("MPIEnvironment") << "MPI initialized without MPI_THREAD_MULTIPLE (not provided)." << std::endl;
        }

        #endif
    }

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    KRATOS_DETAIL("MPIEnvironment") << "MPI initialize called in rank " << rank << "." << std::endl;
}

void MPIEnvironment::Finalize()
{
    if (!MPIEnvironment::IsFinalized())
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        KRATOS_DETAIL("MPIEnvironment") << "MPI finalize called in rank " << rank << "." << std::endl;
        MPI_Finalize();
    }

}

bool MPIEnvironment::IsInitialized()
{
    int mpi_is_initialized;
    MPI_Initialized(&mpi_is_initialized);
    return (bool)mpi_is_initialized;
}

bool MPIEnvironment::IsFinalized()
{
    int mpi_is_finalized;
    MPI_Finalized(&mpi_is_finalized);
    return (bool)mpi_is_finalized;
}

MPI_Comm MPIEnvironment::GetMPICommunicator(const DataCommunicator& rDataCommunicator)
{
    if (rDataCommunicator.IsDistributed())
    {
        const MPIDataCommunicator& r_mpi_data_comm = static_cast<const MPIDataCommunicator&>(rDataCommunicator);
        return r_mpi_data_comm.GetMPICommunicator();
    }
    else
    {
        return MPI_COMM_SELF;
    }
}

void MPIEnvironment::InitializeKratosParallelEnvironment()
{
    // Define the World DataCommunicator as a wrapper for MPI_COMM_WORLD and make it the default.
    ParallelEnvironment& parallel_environment = ParallelEnvironment::GetInstance();
    parallel_environment.RegisterDataCommunicator("World", MPIDataCommunicator(MPI_COMM_WORLD), ParallelEnvironment::MakeDefault);
}

}