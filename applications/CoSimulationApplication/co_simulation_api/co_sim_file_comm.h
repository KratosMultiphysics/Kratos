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
        const SettingsType default_settings = {
            {"communication_folder_name_suffix", ""},
            {"use_folder_for_communication" , "0"}
        };
        Tools::AddMissingSettings(default_settings, CoSimComm::mrSettings);

        mCommFolderSuffix = CoSimComm::mrSettings.at("communication_folder_name_suffix");
        mCommInFolder = (CoSimComm::mrSettings.at("use_folder_for_communication") == "1");
    }

private:

    std::string mCommFolderSuffix = "";
    bool mCommInFolder = false;

    bool ConnectDetail() override
    {
        return true; // nothing needed here for file-based communication
    }

    bool DisconnectDetail() override
    {
        return true; // nothing needed here for file-based communication
    }

    bool ImportDetail(DataContainers::Mesh& rDataContainer, const std::string& rIdentifier) override
    {
        return true;
    }

    bool ExportDetail(const DataContainers::Mesh& rDataContainer, const std::string& rIdentifier) override
    {
        return true;
    }

    bool ImportDetail(DataContainers::Data& rDataContainer, const std::string& rIdentifier) override
    {
        return true;
    }

    bool ExportDetail(const DataContainers::Data& rDataContainer, const std::string& rIdentifier) override
    {
        return true;
    }

};

} // namespace CoSim

#endif /* KRATOS_CO_SIM_FILE_COMM_H_INCLUDED */