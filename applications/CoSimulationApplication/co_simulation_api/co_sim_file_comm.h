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
    explicit FileComm(const std::string& rName, SettingsType& rSettings)
        : CoSimComm(rName, rSettings)
    {
        // throw std::runtime_error("Files Communication is not implemented yet");
    }

private:

    bool ConnectDetail() override
    {
        std::cout << "FileComm; ConnectDetail" << std::endl;
        return true; // nothing needed here for file-based communication
    }

    bool DisconnectDetail() override
    {
        std::cout << "FileComm; DisconnectDetail" << std::endl;

        return true; // nothing needed here for file-based communication
    }

    bool ExportDetail(const DataContainers::Mesh& rDataContainer, const std::string& rIdentifier) override
    {
        return true;
    }

};

} // namespace CoSim

#endif /* KRATOS_CO_SIM_FILE_COMM_H_INCLUDED */