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

#ifndef KRATOS_CO_SIM_API_IMPL_H_INCLUDED
#define KRATOS_CO_SIM_API_IMPL_H_INCLUDED

// Optional includes
#ifdef KRATOS_CO_SIM_API_ENABLE_SOCKETS
#include "co_sim_sockets_comm.h"
#endif /* KRATOS_CO_SIM_API_ENABLE_SOCKETS */


#ifdef KRATOS_CO_SIM_API_ENABLE_MPI
#include "co_sim_mpi_comm.h"
#endif /* KRATOS_CO_SIM_API_ENABLE_MPI */

// System includes
#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <stdexcept>

// Project includes
#include "co_sim_file_comm.h"


namespace CoSim {

CoSimAPI::CoSimAPI(SettingsType& rSettings)
{
    Initialize(rSettings);
}

CoSimAPI::CoSimAPI(const std::string& rSettingsFileName)
{
    // read file
    // Initialize(rSettings);
}

CoSimAPI::~CoSimAPI()
{
    if (!mpComm->Disconnect()) {
        std::cout << "Warning: Disconnect was not successful!" << std::endl;
    }
}


template<class DataContainer>
bool CoSimAPI::SendData(const DataContainer& rContainer, const std::string& rIdentifier)
{
    return mpComm->SendData(rContainer, rIdentifier);
}

template<class DataContainer>
bool CoSimAPI::RecvData(DataContainer& rContainer, const std::string& rIdentifier)
{
    return mpComm->RecvData(rContainer, rIdentifier);
}


void CoSimAPI::Initialize(SettingsType& rSettings)
{
    const std::string comm_format(rSettings.at("communication_format")); // TODO check if specified, if not set to file

    if (comm_format == "file") {
        mpComm = std::unique_ptr<CoSimComm>(new FileComm(rSettings));
    } else if (comm_format == "sockets") {
#ifdef KRATOS_CO_SIM_API_ENABLE_SOCKETS
        mpComm = std::unique_ptr<CoSimComm>(new SocketsComm(rSettings));
#else
        throw std::runtime_error("Support for Sockets was not compiled!");
#endif /* KRATOS_CO_SIM_API_ENABLE_SOCKETS */
    } else if (comm_format == "mpi") {
#ifdef KRATOS_CO_SIM_API_ENABLE_MPI
        mpComm = std::unique_ptr<CoSimComm>(new MPIComm(rSettings));
#else
        throw std::runtime_error("Support for MPI was not compiled!");
#endif /* KRATOS_CO_SIM_API_ENABLE_MPI */
    } else {
        throw std::runtime_error("comm format not available!"); // TODO improve message
    }

    if (!mpComm->Connect()) {
        throw std::runtime_error("Connection was not successful!");
    }
}

} // namespace CoSim

#endif /* KRATOS_CO_SIM_API_IMPL_H_INCLUDED */
