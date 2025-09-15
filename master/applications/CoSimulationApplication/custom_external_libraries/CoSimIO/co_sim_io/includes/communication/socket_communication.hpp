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

#ifndef CO_SIM_IO_SOCKET_COMMUNICATION_INCLUDED
#define CO_SIM_IO_SOCKET_COMMUNICATION_INCLUDED

// System includes

// Project includes
#include "includes/communication/base_socket_communication.hpp"

namespace CoSimIO {
namespace Internals {

class CO_SIM_IO_API SocketCommunication : public BaseSocketCommunication<asio::ip::tcp::socket>
{
public:
    using BaseType = BaseSocketCommunication<asio::ip::tcp::socket>;

    SocketCommunication(
        const Info& I_Settings,
        std::shared_ptr<DataCommunicator> I_DataComm);

    Info ConnectDetail(const Info& I_Info) override;

private:
    std::shared_ptr<asio::ip::tcp::acceptor> mpAsioAcceptor; // probably sufficient to have local in function
    unsigned short mPortNumber=0;
    std::string mIpAddress;
    std::string mSerializedConnectionInfo;

    std::string GetCommunicationName() const override {return "socket";}

    void PrepareConnection(const Info& I_Info) override;

    Info GetCommunicationSettings() const override;

    void GetConnectionInformation();
};

} // namespace Internals
} // namespace CoSimIO

#endif // CO_SIM_IO_SOCKET_COMMUNICATION_INCLUDED
