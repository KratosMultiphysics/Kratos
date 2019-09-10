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

#ifndef KRATOS_CO_SIM_FILE_COMM_H_INCLUDED
#define KRATOS_CO_SIM_FILE_COMM_H_INCLUDED

// System includes

// Project includes
#include "co_sim_comm.h"

namespace CoSim {

class FileComm : public CoSimComm
{
public:
    FileComm(SettingsType& rSettings);

    bool SendData(const DataContainers::Mesh& rContainer, const std::string& rIdentifier) override
    {
        return true;
    }

};

} // namespace CoSim

#endif /* KRATOS_CO_SIM_FILE_COMM_H_INCLUDED */