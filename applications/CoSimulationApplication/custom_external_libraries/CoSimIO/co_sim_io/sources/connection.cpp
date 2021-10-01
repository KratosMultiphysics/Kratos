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

// System includes
#include <iostream>

// Project includes
#include "includes/define.hpp"
#include "includes/connection.hpp"
#include "includes/communication/file_communication.hpp"
#include "includes/communication/pipe_communication.hpp"
#include "includes/communication/sockets_communication.hpp"

namespace CoSimIO {
namespace Internals {

Connection::Connection(
    const Info& I_Settings,
    std::shared_ptr<DataCommunicator> I_DataComm)
    : mpDatacomm(I_DataComm)
{
    Initialize(I_Settings);
}

Info Connection::Connect(const Info& I_Info)
{
    Info info = mpComm->Connect(I_Info);
    info.Set<int>("connection_status", ConnectionStatus::Connected);
    return info;
}

Info Connection::Disconnect(const Info& I_Info)
{
    Info info = mpComm->Disconnect(I_Info);
    info.Set<int>("connection_status", ConnectionStatus::Disconnected);
    return info;
}

Info Connection::Register(
    const std::string& rFunctionName,
    FunctionPointerType FunctionPointer)
{
    CO_SIM_IO_INFO("CoSimIO") << "Registering function for: " << rFunctionName << std::endl;

    CheckIfNameIsValid(rFunctionName);

    CO_SIM_IO_ERROR_IF((mRegisteredFunctions.count(rFunctionName)>0)) << "A function was already registered for " << rFunctionName << "!" << std::endl;

    mRegisteredFunctions[rFunctionName] = FunctionPointer;
    return Info(); // TODO use this
}

Info Connection::Run(const Info& I_Info)
{
    CoSimIO::Info ctrl_info;
    ctrl_info.Set("identifier", "run_control");

    while (true) {
        CoSimIO::Info info = ImportInfo(ctrl_info);
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

void Connection::Initialize(const Info& I_Settings)
{
    const std::string comm_format = I_Settings.Get<std::string>("communication_format", "file"); // default is file-communication

    CO_SIM_IO_INFO_IF("CoSimIO", mpDatacomm->Rank()==0) << "CoSimIO from \"" << I_Settings.Get<std::string>("my_name") << "\" to \"" << I_Settings.Get<std::string>("connect_to") << "\" uses communication format: " << comm_format << std::endl;

    if (comm_format == "file") {
        mpComm = CoSimIO::make_unique<FileCommunication>(I_Settings, mpDatacomm);
    } else if (comm_format == "pipe") {
        mpComm = CoSimIO::make_unique<PipeCommunication>(I_Settings, mpDatacomm);
    } else if (comm_format == "sockets") {
        mpComm = CoSimIO::make_unique<SocketsCommunication>(I_Settings, mpDatacomm);
    } else {
        CO_SIM_IO_ERROR << "Unsupported communication format: " << comm_format << std::endl;
    }
}

void Connection::CheckIfNameIsValid(const std::string& rName) const
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

} // namespace Internals
} // namespace CoSimIO
