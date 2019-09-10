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


#ifdef KRATOS_CO_SIM_API_ENABLE_SOCKETS
#include "sockets_utils.h"
#endif /* KRATOS_CO_SIM_API_ENABLE_SOCKETS */


#ifdef KRATOS_CO_SIM_API_ENABLE_MPI
#include "mpi_utils.h"
#endif /* KRATOS_CO_SIM_API_ENABLE_MPI */

#include <vector>
#include <string>
#include <unordered_map>
#include <stdexcept>

CoSimAPI::CoSimAPI(CoSimAPIConfigType& rConfiguration)
{
    ParseConfiguration(rConfiguration);
}

CoSimAPI::CoSimAPI(const std::string& rConfigurationFileName)
{
    // read file
    // ParseConfiguration(rConfiguration);
}



bool CoSimAPI::SendData(const std::vector<double>& rData, const std::string rIdentifier)
{
    return true;
}

bool CoSimAPI::RecvData(std::vector<double>& rData, const std::string rIdentifier)
{
    return true;
}

void CoSimAPI::ParseConfiguration(CoSimAPIConfigType& rConfiguration)
{
    AssignAndValidateDefaults(rConfiguration);

    mEchoLevel = std::stoi(rConfiguration.at("echo_level"));
    mCommunicationFormat = rConfiguration.at("communication_format");

    if (mCommunicationFormat == "sockets") {
#ifndef KRATOS_CO_SIM_API_ENABLE_SOCKETS
        throw std::runtime_error("Support for Sockets was not compiled!");
#endif /* KRATOS_CO_SIM_API_ENABLE_SOCKETS */
    }

    if (mCommunicationFormat == "mpi") {
#ifndef KRATOS_CO_SIM_API_ENABLE_MPI
        throw std::runtime_error("Support for MPI was not compiled!");
#endif /* KRATOS_CO_SIM_API_ENABLE_MPI */
    }

}

void CoSimAPI::AssignAndValidateDefaults(CoSimAPIConfigType& rConfiguration)
{

}

#endif /* KRATOS_CO_SIM_API_IMPL_H_INCLUDED */
