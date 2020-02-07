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

#include "mpi.h"
#include "mpi/includes/mpi_manager.h"

namespace Kratos
{
    
MPIManager::~MPIManager()
{
    if (!IsFinalized())
    {
        MPI_Finalize();
    }
}

MPIManager::Pointer MPIManager::Create()
{
    return MPIManager::Pointer(new MPIManager());
}

bool MPIManager::IsInitialized() const
{
    int mpi_is_initialized;
    MPI_Initialized(&mpi_is_initialized);
    return (bool)mpi_is_initialized;
}

bool MPIManager::IsFinalized() const
{
    int mpi_is_finalized;
    MPI_Finalized(&mpi_is_finalized);
    return (bool)mpi_is_finalized;
}

MPIManager::MPIManager()
{
    if (!IsInitialized())
    {
        // Note from openmpi docs: https://www.open-mpi.org/doc/v2.0/man3/MPI_Init.3.php
        // "Open MPI accepts the C/C++ argc and argv arguments to main, but neither modifies, interprets, nor distributes them"
        // I am passing empty arguments for simplicity.
        int argc = 0;
        char** argv = nullptr;

        #if MPI_VERSION < 2
        MPI_Init(&argc, &argv);

        KRATOS_DETAIL("MPIManager") << "MPI version older than 2, no support for MPI_THREAD_MULTIPLE is provided." << std::endl;

        #else
        int provided;
        MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

        if (provided < MPI_THREAD_MULTIPLE)
        {
            KRATOS_DETAIL("MPIManager") << "MPI initialized without MPI_THREAD_MULTIPLE (not provided)." << std::endl;
        }

        #endif
    }
}

}