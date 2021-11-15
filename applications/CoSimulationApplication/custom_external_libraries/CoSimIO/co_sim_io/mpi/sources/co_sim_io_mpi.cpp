//     ______     _____ _           ________
//    / ____/___ / ___/(_)___ ___  /  _/ __ |
//   / /   / __ \\__ \/ / __ `__ \ / // / / /
//  / /___/ /_/ /__/ / / / / / / // // /_/ /
//  \____/\____/____/_/_/ /_/ /_/___/\____/
//  Kratos CoSimulationApplication
//
//  License:         BSD License, see license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

/*
This file contains the implementation of the functions defined in "co_sim_io_mpi.hpp"
*/

// System includes

// Project includes
#include "co_sim_io_mpi.hpp"
#include "includes/connect_impl.hpp"
#include "mpi/includes/mpi_data_communicator.hpp"

namespace CoSimIO {

Info ConnectMPI(
    const Info& I_Settings,
    MPI_Comm ThisMPIComm)
{
    // make sure MPI was already initialized by the host
    int flag_initialized;
    MPI_Initialized(&flag_initialized);
    CO_SIM_IO_ERROR_IF_NOT(flag_initialized) << "MPI must be initialized before calling \"ConnectMPI\"!" << std::endl;

    // make sure a valid communicator is passed
    CO_SIM_IO_ERROR_IF(ThisMPIComm == MPI_COMM_SELF) << "Passing \"MPI_COMM_SELF\" is not allowed!" << std::endl;
    CO_SIM_IO_ERROR_IF(ThisMPIComm == MPI_COMM_NULL) << "Passing \"MPI_COMM_NULL\" is not allowed!" << std::endl;

    return Internals::ConnectImpl(I_Settings, std::make_shared<CoSimIO::Internals::MPIDataCommunicator>(ThisMPIComm));
}

} // namespace CoSimIO
