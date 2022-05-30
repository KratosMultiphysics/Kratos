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

// System includes
#include <unordered_map>
#include <memory>
#include <string>
#include <functional>
#include <utility>

// Project includes
#include "includes/info.hpp"
#include "includes/data_communicator.hpp"
#include "includes/communication/communication.hpp"

namespace CoSimIO {
namespace Internals {

class CO_SIM_IO_API Connection
{

public:

    using FunctionPointerType = std::function<Info(const Info&)>;

    Connection(
        const Info& I_Settings,
        std::shared_ptr<DataCommunicator> I_DataComm);

    Info Connect(const Info& I_Info);

    Info Disconnect(const Info& I_Info);

    Info Register(
        const std::string& rFunctionName,
        FunctionPointerType FunctionPointer);

    Info Run(const Info& I_Info);


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

    std::shared_ptr<DataCommunicator> mpDatacomm;

    std::unordered_map<std::string, FunctionPointerType> mRegisteredFunctions;

    void Initialize(const Info& I_Settings);

    void CheckIfNameIsValid(const std::string& rName) const;

}; // class Connection

} // namespace Internals
} // namespace CoSimIO

#endif // CO_SIM_IO_CONNECTION_INCLUDED
