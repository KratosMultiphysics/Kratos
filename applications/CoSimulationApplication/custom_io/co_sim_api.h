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

#ifndef KRATOS_CO_SIM_API_H_INCLUDED
#define KRATOS_CO_SIM_API_H_INCLUDED

#define KRATOS_CO_SIM_API_VERSION 1 // TODO probably this has to be in a different file ...

/*
This file defines the API of Kratos-CoSimulation for the exchange of data
with external solvers.
By default the communication is done through files,
support for sockets and MPI can optionally be enabled
*/

// #define KRATOS_CO_SIM_API_ENABLE_SOCKETS // uncomment for Sockets support
// #define KRATOS_CO_SIM_API_ENABLE_MPI // uncomment for MPI support

// System includes
#include <string>
#include <memory>

// Project includes
#include "co_sim_comm.h"

namespace CoSim {

class CoSimAPI
{

public:
    typedef CoSimComm::SettingsType SettingsType;

    // Constructor establishes connection
    CoSimAPI(SettingsType& rSettings);
    CoSimAPI(const std::string& rSettingsFileName);

    // Destructor performs disconnect
    ~CoSimAPI();

    template<class DataContainer>
    bool SendData(const DataContainer& rContainer, const std::string& rIdentifier);

    template<class DataContainer>
    bool RecvData(DataContainer& rContainer, const std::string& rIdentifier);

private:
    std::unique_ptr<CoSimComm> mpComm;

    void Initialize(SettingsType& rSettings);

}; // class CoSimAPI

} // namespace CoSim

#include "co_sim_api_impl.h"

#endif /* KRATOS_CO_SIM_API_H_INCLUDED */
