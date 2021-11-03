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

#ifndef CO_SIM_IO_MPI_INCLUDED
#define CO_SIM_IO_MPI_INCLUDED

/*
This file contains the interface of the MPI extension for the CoSimIO
*/

// System includes
#include "mpi.h"

// Project includes
#include "co_sim_io.hpp"

namespace CoSimIO {

Info ConnectMPI(
    const Info& I_Settings,
    MPI_Comm ThisMPIComm);

} // namespace CoSimIO

#endif // CO_SIM_IO_MPI_INCLUDED
