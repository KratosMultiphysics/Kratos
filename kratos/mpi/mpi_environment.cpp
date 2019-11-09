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
#ifdef _OPENMP
#include "omp.h"
#endif

// Project includes
#include "input_output/logger.h"
#include "includes/parallel_environment.h"

// Module includes
#include "mpi_environment.h"
#include "mpi/includes/mpi_data_communicator.h"

namespace Kratos {

MPIEnvironment& MPIEnvironment::Instance()
{
    // Using double-checked locking to ensure thread safety in the first creation of the singleton.
    if (mpInstance == nullptr)
    {
        #ifdef _OPENMP
        #pragma omp critical
        if (mpInstance == nullptr)
        {
        #endif
            KRATOS_ERROR_IF(mDestroyed) << "Accessing MPIEnvironment after its destruction" << std::endl;
            Create();
        #ifdef _OPENMP
        }
        #endif
    }

    return *mpInstance;
}

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

        KRATOS_DETAIL("MPIEnvironment") << "MPI version older than 2, no support for MPI_THREAD_MULTIPLE is provided." << std::endl;

        #else
        int provided;
        MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

        if (provided < MPI_THREAD_MULTIPLE)
        {
            KRATOS_DETAIL("MPIEnvironment") << "MPI initialized without MPI_THREAD_MULTIPLE (not provided)." << std::endl;
        }

        #endif
    }

    KRATOS_DETAIL_ALL_RANKS("MPIEnvironment") << "MPI initialize called." << std::endl;
}

void MPIEnvironment::Finalize()
{
    if (!MPIEnvironment::IsFinalized())
    {
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

MPIEnvironment::MPIEnvironment() {}

MPIEnvironment::~MPIEnvironment()
{
    auto& r_instance = Instance();
    r_instance.Finalize();
    mpInstance = nullptr;
    mDestroyed = true;
}

void MPIEnvironment::Create()
{
    static MPIEnvironment mpi_environment;
    mpInstance = &mpi_environment;
}

MPIEnvironment* MPIEnvironment::mpInstance = nullptr;
bool MPIEnvironment::mDestroyed = false;

}