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

#ifndef CO_SIM_IO_BASE_SOCKET_COMMUNICATION_INCLUDED
#define CO_SIM_IO_BASE_SOCKET_COMMUNICATION_INCLUDED

// System includes
#include <thread>

// Project includes
#include "communication.hpp"

// External includes
#define ASIO_NO_DEPRECATED // disabling deprecated features/interfaces
#define ASIO_STANDALONE // independent of boost
#ifndef _WIN32_WINNT
    #define _WIN32_WINNT 0x0601 // see "https://github.com/chriskohlhoff/asio/issues/596"
#endif
#include "asio.hpp"

namespace CoSimIO {
namespace Internals {

template<class TSocketType>
class CO_SIM_IO_API BaseSocketCommunication : public Communication
{
public:
    BaseSocketCommunication(
        const Info& I_Settings,
        std::shared_ptr<DataCommunicator> I_DataComm)
        : Communication(I_Settings, I_DataComm) {}

    ~BaseSocketCommunication() override;

    Info ConnectDetail(const Info& I_Info) override;

    Info DisconnectDetail(const Info& I_Info) override;

protected:
    std::shared_ptr<TSocketType> mpAsioSocket;
    asio::io_context mAsioContext;
    std::thread mContextThread;

    double SendString(
        const Info& I_Info,
        const std::string& rData) override;

    double ReceiveString(
        const Info& I_Info,
        std::string& rData) override;

    double SendDataContainer(
        const Info& I_Info,
        const Internals::DataContainer<double>& rData) override;

    double ReceiveDataContainer(
        const Info& I_Info,
        Internals::DataContainer<double>& rData) override;

    void SendSize(const std::uint64_t Size);

    std::uint64_t ReceiveSize();
};

} // namespace Internals
} // namespace CoSimIO

#endif // CO_SIM_IO_BASE_SOCKET_COMMUNICATION_INCLUDED
