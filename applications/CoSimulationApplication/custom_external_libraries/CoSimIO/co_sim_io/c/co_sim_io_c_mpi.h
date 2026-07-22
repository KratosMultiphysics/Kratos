/*   ______     _____ _           ________
    / ____/___ / ___/(_)___ ___  /  _/ __ |
   / /   / __ \\__ \/ / __ `__ \ / // / / /
  / /___/ /_/ /__/ / / / / / / // // /_/ /
  \____/\____/____/_/_/ /_/ /_/___/\____/
  Kratos CoSimulationApplication

  License:         BSD License, see license.txt

  Main authors:    Philipp Bucher (https://github.com/philbucher)
*/

#ifndef CO_SIM_IO_C_MPI_INCLUDED
#define CO_SIM_IO_C_MPI_INCLUDED

/* C MPI-Interface for CoSimulation
   see "co_sim_io_mpi.hpp"
*/

#include "mpi.h"
#include "co_sim_io_c.h"

#ifdef __cplusplus
extern "C" { /* Define extern C if C++ compiler is used */
#endif

CO_SIM_IO_NODISCARD CoSimIO_Info CoSimIO_ConnectMPI(
    const CoSimIO_Info I_Settings,
    MPI_Comm ThisMPIComm);

#ifdef __cplusplus
}
#endif

#endif /* CO_SIM_IO_C_MPI_INCLUDED */
