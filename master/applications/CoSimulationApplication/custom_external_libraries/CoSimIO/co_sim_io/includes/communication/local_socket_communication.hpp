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

#ifndef CO_SIM_IO_LOCAL_SOCKET_COMMUNICATION_INCLUDED
#define CO_SIM_IO_LOCAL_SOCKET_COMMUNICATION_INCLUDED

// System includes

// Project includes
#include "includes/communication/base_socket_communication.hpp"

namespace CoSimIO {
namespace Internals {

class CO_SIM_IO_API LocalSocketCommunication : public BaseSocketCommunication<asio::local::stream_protocol::socket>
{
public:
    using BaseType = BaseSocketCommunication<asio::local::stream_protocol::socket>;

    LocalSocketCommunication(
        const Info& I_Settings,
        std::shared_ptr<DataCommunicator> I_DataComm)
        : BaseType(I_Settings, I_DataComm) {}

    Info ConnectDetail(const Info& I_Info) override;

private:
    std::shared_ptr<asio::local::stream_protocol::acceptor> mpAsioAcceptor; // probably sufficient to have local in function

    std::string GetCommunicationName() const override {return "local_socket";}
};

} // namespace Internals
} // namespace CoSimIO

#endif // CO_SIM_IO_LOCAL_SOCKET_COMMUNICATION_INCLUDED
