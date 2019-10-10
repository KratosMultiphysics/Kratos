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
#include <unordered_map>
#include <string>
#include <stdexcept>

// Project includes
#include "co_sim_data_containers.h"

namespace CoSim {

#define CO_SIM_COMM_REGISTER_DATA_CONTAINER_TYPE(DataContainerType)                              \
    virtual bool Import(DataContainerType& rDataContainer, const std::string& rIdentifier)       \
        { throw std::runtime_error("Type of data not yet supported"); }                          \
    virtual bool Export(const DataContainerType& rDataContainer, const std::string& rIdentifier) \
        { throw std::runtime_error("Type of data not yet supported"); }                          \

class CoSimComm
{
public:
    typedef std::unordered_map<std::string, std::string> SettingsType;

    virtual bool Connect()    { return true; };
    virtual bool Disconnect() { return true; };

    CO_SIM_COMM_REGISTER_DATA_CONTAINER_TYPE(std::vector<double>);
    CO_SIM_COMM_REGISTER_DATA_CONTAINER_TYPE(DataContainers::Mesh);
    CO_SIM_COMM_REGISTER_DATA_CONTAINER_TYPE(DataContainers::ConvergenceSignal);

private:
    void AssignAndValidateDefaults(const SettingsType& rDefaultSettings, SettingsType& rSettings)
    {
        // ...
    }

};

#undef CO_SIM_COMM_REGISTER_DATA_CONTAINER_TYPE // this macro should only be used for this class

} // namespace CoSim

#endif /* KRATOS_CO_SIM_COMM_H_INCLUDED */
