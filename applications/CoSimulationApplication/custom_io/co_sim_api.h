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

/*
This file defines the API of Kratos-CoSimulation for the exchange of data
with external solvers.
By default the communication is done through files,
support for sockets and MPI can optionally be enabled
*/

// #define KRATOS_CO_SIM_API_ENABLE_SOCKETS // uncomment for Sockets support
// #define KRATOS_CO_SIM_API_ENABLE_MPI // uncomment for MPI support

#include <vector>
#include <string>
#include <unordered_map>

typedef std::unordered_map<std::string, std::string> CoSimAPIConfigType;

// // stores the coordinates in an array [x1,y1,z1,x2,y2,z2,x3,y3,z3,...]
// struct Coordinates
// {

// };

// struct Connectivities
// {

// };

// struct Ids
// {

// };

// struct ConvergenceSignal
// {

// };

class CoSimAPI
{

public:
    CoSimAPI(CoSimAPIConfigType& rConfiguration);
    CoSimAPI(const std::string& rConfigurationFileName);


    bool SendData(const std::vector<double>& rData, const std::string rIdentifier);

    // template<class DataContainer>
    // bool SendData(const DataContainer& rContainer,  const std::string rIdentifier);


    bool RecvData(std::vector<double>& rData,       const std::string rIdentifier);

    // template<class DataContainer>
    // bool RecvData(const DataContainer& rContainer,  const std::string rIdentifier);

private:
    std::string mCommunicationFormat;
    int mEchoLevel;

    CoSimAPIConfigType mDefaults {
        { "communication_format", "file" },
        { "echo_level", "0" }
    };

    void AssignAndValidateDefaults(CoSimAPIConfigType& rConfiguration);

    void ParseConfiguration(CoSimAPIConfigType& rConfiguration);

}; // class CoSimAPI

#include "co_sim_api_impl.h"

#endif /* KRATOS_CO_SIM_API_H_INCLUDED */
