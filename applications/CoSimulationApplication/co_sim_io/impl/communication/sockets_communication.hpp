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

#ifndef CO_SIM_IO_SOCKETS_COMMUNICATION_H_INCLUDED
#define CO_SIM_IO_SOCKETS_COMMUNICATION_H_INCLUDED

// System includes

// External includes
#include <boost/asio.hpp>

// Project includes
#include "communication.hpp"

namespace CoSimIO {
namespace Internals {

class SocketsCommunication : public Communication
{
public:
    explicit SocketsCommunication(const std::string& rName, SettingsType& rSettings, const bool IsConnectionMaster)
        : Communication(rName, rSettings, IsConnectionMaster)
    {
       static_assert(false,"Sockets Communication is not implemented yet");
    }
};

} // namespace Internals
} // namespace CoSimIO

#endif // CO_SIM_IO_SOCKETS_COMMUNICATION_H_INCLUDED
