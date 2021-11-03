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

// External includes
#define ASIO_NO_DEPRECATED // disabling deprecated features/interfaces
#define ASIO_STANDALONE // independent of boost
#ifndef _WIN32_WINNT
    #define _WIN32_WINNT 0x0601 // see "https://github.com/chriskohlhoff/asio/issues/596"
#endif
#include "asio.hpp"

// Project includes
#include "includes/communication/sockets_communication.hpp"

namespace CoSimIO {
namespace Internals {

SocketsCommunication::SocketsCommunication(
    const Info& I_Settings,
    std::shared_ptr<DataCommunicator> I_DataComm)
    : Communication(I_Settings, I_DataComm)
{
    CO_SIM_IO_ERROR << "Sockets Communication is not implemented yet" << std::endl;
}

Info SocketsCommunication::ConnectDetail(const Info& I_Info)
{
    return Info();
}

Info SocketsCommunication::DisconnectDetail(const Info& I_Info)
{
    return Info();
}

} // namespace Internals
} // namespace CoSimIO
