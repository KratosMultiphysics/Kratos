//     ______     _____ _           ________
//    / ____/___ / ___/(_)___ ___  /  _/ __ |
//   / /   / __ \\__ \/ / __ `__ \ / // / / /
//  / /___/ /_/ /__/ / / / / / / // // /_/ /
//  \____/\____/____/_/_/ /_/ /_/___/\____/
//  Kratos CoSimulationApplication
//
//  License:         BSD License, see license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

#ifndef CO_SIM_IO_CONNECTION_INCLUDED
#define CO_SIM_IO_CONNECTION_INCLUDED

// Optional includes
#ifdef CO_SIM_IO_USING_SOCKETS
#include "communication/sockets_communication.hpp"
#endif // CO_SIM_IO_USING_SOCKETS


#ifdef CO_SIM_IO_USING_MPI
#include "communication/mpi_communication.hpp"
#endif // CO_SIM_IO_USING_MPI

// System includes
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <utility>
#include <functional>

// Project includes
#include "communication/file_communication.hpp"

namespace CoSimIO {
namespace Internals {

class Connection
{

public:

    using FunctionPointerType = std::function<Info(const Info&)>;

    explicit Connection(const Info& I_Settings)
    {
        Initialize(I_Settings);
    }

    Info Connect(const Info& I_Info)
    {
        Info info = mpComm->Connect(I_Info);
        info.Set<int>("connection_status", ConnectionStatus::Connected);
        return info;
    }

    Info Disconnect(const Info& I_Info)
    {
        Info info = mpComm->Disconnect(I_Info);
        info.Set<int>("connection_status", ConnectionStatus::Disconnected);
        return info;
    }

    Info Register(
        const std::string& rFunctionName,
        FunctionPointerType FunctionPointer)
    {
        CO_SIM_IO_INFO("CoSimIO") << "Registering function for: " << rFunctionName << std::endl;

        CheckIfNameIsValid(rFunctionName);

        CO_SIM_IO_ERROR_IF((mRegisteredFunctions.count(rFunctionName)>0)) << "A function was already registered for " << rFunctionName << "!" << std::endl;

        mRegisteredFunctions[rFunctionName] = FunctionPointer;
        return Info(); // TODO use this
    }

    Info Run(const Info& I_Info)
    {
        CoSimIO::Info ctrl_info;
        ctrl_info.Set("identifier", "run_control");

        while(true) {
            Info info = ImportInfo(ctrl_info);
            const std::string control_signal = info.Get<std::string>("control_signal");
            CheckIfNameIsValid(control_signal);
            if (control_signal == "exit") {
                break;
            } else {
                auto it_fct = mRegisteredFunctions.find(control_signal);
                if (it_fct == mRegisteredFunctions.end()) {
                    std::stringstream err_msg;
                    err_msg << "Nothing was registered for \"" << control_signal << "\"!\nOnly the following names are currently registered:";
                    for (const auto& reg : mRegisteredFunctions) {
                        err_msg << "\n    " << reg.first;
                    }
                    err_msg << "\n    end" << std::endl;
                    CO_SIM_IO_ERROR << err_msg.str();
                }
                it_fct->second(info.Get<Info>("settings", Info{})); // pass settings if specified

            }
        }
        return Info(); // TODO use this
    }


    template<class... Args>
    Info ImportInfo(Args&&... args)
    {
        return mpComm->ImportInfo(std::forward<Args>(args)...);
    }

    template<class... Args>
    Info ExportInfo(Args&&... args)
    {
        return mpComm->ExportInfo(std::forward<Args>(args)...);
    }
    template<class... Args>
    Info ImportData(Args&&... args)
    {
        return mpComm->ImportData(std::forward<Args>(args)...);
    }

    template<class... Args>
    Info ExportData(Args&&... args)
    {
        return mpComm->ExportData(std::forward<Args>(args)...);
    }

    template<class... Args>
    Info ImportMesh(Args&&... args)
    {
        return mpComm->ImportMesh(std::forward<Args>(args)...);
    }

    template<class... Args>
    Info ExportMesh(Args&&... args)
    {
        return mpComm->ExportMesh(std::forward<Args>(args)...);
    }

private:
    std::unique_ptr<Communication> mpComm; // handles communication (File, Sockets, MPI, ...)

    std::unordered_map<std::string, FunctionPointerType> mRegisteredFunctions;

    void Initialize(const Info& I_Settings)
    {
        const std::string comm_format = I_Settings.Get<std::string>("communication_format", "file"); // default is file-communication

        CO_SIM_IO_INFO("CoSimIO") << "CoSimIO from \"" << I_Settings.Get<std::string>("my_name") << "\" to \"" << I_Settings.Get<std::string>("connect_to") << "\" uses communication format: " << comm_format << std::endl;

        if (comm_format == "file") {
            mpComm = CoSimIO::make_unique<FileCommunication>(I_Settings);
        } else if (comm_format == "sockets") {
            #ifdef CO_SIM_IO_USING_SOCKETS
            mpComm = CoSimIO::make_unique<SocketsCommunication>(I_Settings);
            #else
            CO_SIM_IO_ERROR << "Support for Sockets was not compiled!" << std::endl;
            #endif // CO_SIM_IO_USING_SOCKETS
        } else if (comm_format == "mpi") {
            #ifdef CO_SIM_IO_USING_MPI
            mpComm = CoSimIO::make_unique<MPICommunication>(I_Settings);
            #else
            CO_SIM_IO_ERROR << "Support for MPI was not compiled!" << std::endl;
            #endif // CO_SIM_IO_USING_MPI
        } else {
            CO_SIM_IO_ERROR << "Unsupported communication format: " << comm_format << std::endl;
        }
    }

    void CheckIfNameIsValid(const std::string& rName) const
    {
        // could use set but that would require another include just for this
        const static std::vector<std::string> allowed_names {
            "AdvanceInTime",
            "InitializeSolutionStep",
            "Predict",
            "SolveSolutionStep",
            "FinalizeSolutionStep",
            "OutputSolutionStep",
            "ImportMesh",
            "ExportMesh",
            "ImportData",
            "ExportData",
            "exit"
        };

        if (std::find(allowed_names.begin(), allowed_names.end(), rName) == allowed_names.end()) {
            std::stringstream err_msg;
            err_msg << "The name \"" << rName << "\" is not allowed!\nOnly the following names are allowed:";
            for (const auto& name : allowed_names) {
                err_msg << "\n    " << name;
            }
            err_msg << std::endl;
            CO_SIM_IO_ERROR << err_msg.str();
        }
    }

}; // class Connection

} // namespace Internals
} // namespace CoSimIO

#endif // CO_SIM_IO_CONNECTION_INCLUDED
