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

CoSimIO::CoSimIO(const std::string rName, const std::string& rSettingsFileName)
    : CoSimIO::CoSimIO(rName, ReadSettingsFile(rSettingsFileName)) { } // forwarding constructor call

CoSimIO::CoSimIO(const std::string rName, SettingsType rSettings)
{
    Initialize(rSettings);
}

bool CoSimIO::Connect()
{
    return mpComm->Connect();
}

bool CoSimIO::Disconnect()
{
    return mpComm->Disconnect();
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


void CoSimIO::Initialize(SettingsType& rSettings)
{
    AddMissingSettings(rSettings);

    const std::string comm_format(rSettings.at("communication_format"));
    mEchoLevel = std::stoi(rSettings.at("echo_level"));

    std::cout << "CoSimIO uses the following configuration:";
    std::cout << "\n    Communication Format: " << comm_format;
    std::cout << "\n    Echo Level: " << mEchoLevel << std::endl;

    if (comm_format == "file") {
        mpComm = std::unique_ptr<CoSimComm>(new FileComm(rSettings)); // make_unique is C++14
    } else if (comm_format == "sockets") {
#ifdef KRATOS_CO_SIM_IO_ENABLE_SOCKETS
        mpComm = std::unique_ptr<CoSimComm>(new SocketsComm(rSettings)); // make_unique is C++14
#else
        throw std::runtime_error("Support for Sockets was not compiled!");
#endif /* KRATOS_CO_SIM_IO_ENABLE_SOCKETS */
    } else if (comm_format == "mpi") {
#ifdef KRATOS_CO_SIM_IO_ENABLE_MPI
        mpComm = std::unique_ptr<CoSimComm>(new MPIComm(rSettings)); // make_unique is C++14
#else
        throw std::runtime_error("Support for MPI was not compiled!");
#endif /* KRATOS_CO_SIM_IO_ENABLE_MPI */
    } else {
        std::stringstream err_msg;
        err_msg << "Unsupported communication format: " << comm_format;
        throw std::runtime_error(err_msg.str());
    }
}

void CoSimIO::AddMissingSettings(SettingsType& rSettings)
{
    const SettingsType default_settings = {
		{"communication_format", "file"},
		{"echo_level",           "1"}
    };

    for (const auto& r_setting : default_settings) {
        if (rSettings.count(r_setting.first) == 0) {
            rSettings[r_setting.first] = r_setting.second;
        }
    }
}

CoSimIO::SettingsType CoSimIO::ReadSettingsFile(const std::string& rSettingsFileName)
{
    std::ifstream settings_file(rSettingsFileName);

    if (!settings_file.good()) {
        std::cout << "Input file \"" << rSettingsFileName << "\" could not be read, using default configuration" << std::endl;
        return SettingsType();
    }

    std::string current_line;

    while (std::getline(settings_file, current_line)) {
        // TODO implement this
    }
}

} // namespace CoSim

#endif /* KRATOS_CO_SIM_IO_IMPL_H_INCLUDED */
