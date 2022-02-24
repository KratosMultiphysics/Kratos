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
#include "includes/communication/factory.hpp"

namespace CoSimIO {
namespace Internals {

Connection::Connection(
    const Info& I_Settings,
    std::shared_ptr<DataCommunicator> I_DataComm,
    const CommunicationFactory& rCommFactory)
    : mpDatacomm(I_DataComm)
{
    Initialize(I_Settings, rCommFactory);
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

void Connection::Initialize(
    const Info& I_Settings,
    const CommunicationFactory& rCommFactory)
{
    Info comm_settings(I_Settings);
    if (!comm_settings.Has("communication_format")) {
        // set default communication format if not provided by user
        // default is socket-communication
        comm_settings.Set<std::string>("communication_format", "socket");
    }

    const std::string comm_format = comm_settings.Get<std::string>("communication_format");

    CO_SIM_IO_INFO_IF("CoSimIO", mpDatacomm->Rank()==0) << "CoSimIO from \"" << comm_settings.Get<std::string>("my_name") << "\" to \"" << comm_settings.Get<std::string>("connect_to") << "\" uses communication format: " << comm_format << std::endl;

    mpComm = rCommFactory.Create(comm_settings, mpDatacomm);
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
        CO_SIM_IO_ERROR << err_msg.str() << std::endl;
    }
}

} // namespace Internals
} // namespace CoSimIO
