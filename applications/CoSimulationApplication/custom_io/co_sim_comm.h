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

#ifndef KRATOS_CO_SIM_COMM_H_INCLUDED
#define KRATOS_CO_SIM_COMM_H_INCLUDED

// System includes
#include <stdexcept>

// Project includes
#include "co_sim_data_containers.h"

namespace CoSim {

class CoSimComm
{
public:
    typedef std::unordered_map<std::string, std::string> SettingsType;

    virtual bool Connect()    { return true; };
    virtual bool Disconnect() { return true; };

    virtual bool SendData(const std::vector<double>& rData, const std::string& rIdentifier)
    {
        throw std::runtime_error("Type of data not yet supported");
    }

    virtual bool SendData(const DataContainers::Mesh& rContainer, const std::string& rIdentifier)
    {
        throw std::runtime_error("Type of data not yet supported");
    }

    virtual bool SendData(const DataContainers::ConvergenceSignal& rContainer, const std::string& rIdentifier)
    {
        throw std::runtime_error("Type of data not yet supported");
    }

    virtual bool RecvData(std::vector<double>& rData, const std::string& rIdentifier)
    {
        throw std::runtime_error("Type of data not yet supported");
    }

    virtual bool RecvData(DataContainers::Mesh& rContainer, const std::string& rIdentifier)
    {
        throw std::runtime_error("Type of data not yet supported");
    }

    virtual bool RecvData(DataContainers::ConvergenceSignal& rContainer, const std::string& rIdentifier)
    {
        throw std::runtime_error("Type of data not yet supported");
    }

private:
    void AssignAndValidateDefaults(const SettingsType& rDefaultSettings, SettingsType& rSettings)
    {
        // ...
    }

};

} // namespace CoSim

#endif /* KRATOS_CO_SIM_COMM_H_INCLUDED */
