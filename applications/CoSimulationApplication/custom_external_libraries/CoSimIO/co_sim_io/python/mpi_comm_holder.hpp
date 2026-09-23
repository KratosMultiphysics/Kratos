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

#ifndef CO_SIM_IO_PYHON_MPI_COMM_HOLDER_INCLUDED
#define CO_SIM_IO_PYHON_MPI_COMM_HOLDER_INCLUDED

// External includes
#include "mpi.h"

namespace CoSimIO {

class MPICommHolder
{
public:
    MPICommHolder(MPI_Comm Comm) : mComm(Comm) {}
    virtual ~MPICommHolder() = default;

    MPI_Comm GetMPIComm() const {return mComm;}

private:
    MPI_Comm mComm;
};

} // namespace CoSimIO

#endif // CO_SIM_IO_PYHON_MPI_COMM_HOLDER_INCLUDED
