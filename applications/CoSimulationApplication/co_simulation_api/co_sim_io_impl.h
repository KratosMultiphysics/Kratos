// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:		 BSD License
//                   license: CoSimulationApplication/license.txt
//
//  Main authors:    Philipp Bucher
//

#ifndef KRATOS_CO_SIM_IO_IMPL_H_INCLUDED
#define KRATOS_CO_SIM_IO_IMPL_H_INCLUDED

// Optional includes
#ifdef KRATOS_CO_SIM_IO_ENABLE_SOCKETS
#include "co_sim_sockets_comm.h"
#endif /* KRATOS_CO_SIM_IO_ENABLE_SOCKETS */


#ifdef KRATOS_CO_SIM_IO_ENABLE_MPI
#include "co_sim_mpi_comm.h"
#endif /* KRATOS_CO_SIM_IO_ENABLE_MPI */

// System includes
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>

// Project includes
#include "co_sim_file_comm.h"


namespace CoSim {

inline CoSimIO::CoSimIO(const std::string& rName, const std::string& rSettingsFileName, const bool IsConnectionMaster)
    : CoSimIO::CoSimIO(rName, Tools::ReadSettingsFile(rSettingsFileName), IsConnectionMaster) { } // forwarding constructor call

inline CoSimIO::CoSimIO(const std::string& rName, SettingsType rSettings, const bool IsConnectionMaster)
{
    Initialize(rName, rSettings, IsConnectionMaster);
}

inline bool CoSimIO::Connect()
{
    return mpComm->Connect();
}

inline bool CoSimIO::Disconnect()
{
    return mpComm->Disconnect();
}


inline void CoSimIO::SendControlSignal(const int rSignal, const std::string& rIdentifier)
{
    mpComm->Export(rSignal, rIdentifier);

}
inline int CoSimIO::RecvControlSignal(std::string& rIdentifier)
{
    int signal;
    mpComm->Import(signal, rIdentifier);
    return signal;
}

template<class DataContainer>
bool CoSimIO::Import(DataContainer& rContainer, const std::string& rIdentifier)
{
    return mpComm->Import(rContainer, rIdentifier);
}

template<class DataContainer>
bool CoSimIO::Export(const DataContainer& rContainer, const std::string& rIdentifier)
{
    return mpComm->Export(rContainer, rIdentifier);
}

inline void CoSimIO::Initialize(const std::string& rName, SettingsType& rSettings, const bool IsConnectionMaster)
{
    std::string comm_format("file"); // default is file-communication
    if (rSettings.count("communication_format") != 0) { // communication format has been specified
        comm_format = rSettings.at("communication_format");
    }

    std::cout << "CoSimIO for \"" << rName << "\" uses communication format: " << comm_format << std::endl;

    if (comm_format == "file") {
        mpComm = std::unique_ptr<CoSimComm>(new FileComm(rName, rSettings, IsConnectionMaster)); // make_unique is C++14
    } else if (comm_format == "sockets") {
#ifdef KRATOS_CO_SIM_IO_ENABLE_SOCKETS
        mpComm = std::unique_ptr<CoSimComm>(new SocketsComm(rName, rSettings, IsConnectionMaster)); // make_unique is C++14
#else
        throw std::runtime_error("Support for Sockets was not compiled!");
#endif /* KRATOS_CO_SIM_IO_ENABLE_SOCKETS */
    } else if (comm_format == "mpi") {
#ifdef KRATOS_CO_SIM_IO_ENABLE_MPI
        mpComm = std::unique_ptr<CoSimComm>(new MPIComm(rName, rSettings, IsConnectionMaster)); // make_unique is C++14
#else
        throw std::runtime_error("Support for MPI was not compiled!");
#endif /* KRATOS_CO_SIM_IO_ENABLE_MPI */
    } else {
        std::stringstream err_msg;
        err_msg << "Unsupported communication format: " << comm_format;
        throw std::runtime_error(err_msg.str());
    }
}

} // namespace CoSim

#endif /* KRATOS_CO_SIM_IO_IMPL_H_INCLUDED */
