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

// Project includes
#include "co_sim_io_c_mpi.h"
#include "co_sim_io_mpi.hpp"

namespace {
    // get C Info from C++ Info
    CoSimIO_Info ConvertInfo(CoSimIO::Info I_Info) {
        CoSimIO_Info info;
        info.PtrCppInfo = new CoSimIO::Info(I_Info);
        return info;
    }

    // get C++ Info from C Info
    CoSimIO::Info ConvertInfo(CoSimIO_Info I_Info) {
        return CoSimIO::Info(*(static_cast<CoSimIO::Info*>(I_Info.PtrCppInfo)));
    }
}

CO_SIM_IO_NODISCARD CoSimIO_Info CoSimIO_ConnectMPI(
    const CoSimIO_Info I_Settings,
    MPI_Comm ThisMPIComm)
{
    return ConvertInfo(CoSimIO::ConnectMPI(ConvertInfo(I_Settings), ThisMPIComm));
}
