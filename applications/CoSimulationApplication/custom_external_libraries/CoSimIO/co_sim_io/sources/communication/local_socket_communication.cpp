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
#include <algorithm>

// Project includes
#include "includes/communication/local_socket_communication.hpp"

namespace CoSimIO {
namespace Internals {

Info LocalSocketCommunication::ConnectDetail(const Info& I_Info)
{
    CO_SIM_IO_TRY

    CO_SIM_IO_INFO_IF("CoSimIO", GetDataCommunicator().IsDistributed() && GetDataCommunicator().Rank()==0) << "Warning: Connection was done with MPI, but local-socket based communication works only within the same machine. Communicating between different compute nodes in a distributed memory machine when does not work, it will hang!" << std::endl;

    using asio::local::stream_protocol;
    mpAsioSocket = std::make_shared<stream_protocol::socket>(mAsioContext);

    const std::string bind_file_name = fs::path(GetCommunicationDirectory() / "socket_bind").string();

    if (GetIsPrimaryConnection()) { // this is the server
        std::ofstream bind_file;
        bind_file.open(bind_file_name);
        bind_file.close();
        #ifdef CO_SIM_IO_COMPILED_IN_WINDOWS
        ::_unlink(bind_file_name.c_str()); // Remove previous binding.
        #else
        ::unlink(bind_file_name.c_str()); // Remove previous binding.
        #endif
    }

    SynchronizeAll("local_sock_1");

    stream_protocol::endpoint this_endpoint(bind_file_name.c_str());

    if (GetIsPrimaryConnection()) { // this is the server
        mpAsioAcceptor = std::make_shared<stream_protocol::acceptor>(mAsioContext, this_endpoint);
        SynchronizeAll("local_sock_2");
        mpAsioAcceptor->accept(*mpAsioSocket);
        mpAsioAcceptor->close();
    } else { // this is the client
        SynchronizeAll("local_sock_2");
        mpAsioSocket->connect(this_endpoint);
    }

    return BaseType::ConnectDetail(I_Info);

    CO_SIM_IO_CATCH
}

} // namespace Internals
} // namespace CoSimIO
