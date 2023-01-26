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
#include "includes/communication/socket_communication.hpp"

#ifdef CO_SIM_IO_COMPILED_IN_LINUX
// to detect network interfaces
#include <ifaddrs.h>
#endif

namespace CoSimIO {
namespace Internals {

namespace {

static std::string LOCAL_IP_ADDRESS {"127.0.0.1"}; // local loopback interface

std::unordered_map<std::string, std::string> GetIpv4Addresses()
{
    CO_SIM_IO_TRY

    std::unordered_map<std::string, std::string> networks;

#ifdef CO_SIM_IO_COMPILED_IN_LINUX
    // adapted from: https://www.cyberithub.com/list-network-interfaces/
    struct ifaddrs *addresses;
    CO_SIM_IO_ERROR_IF(getifaddrs(&addresses) == -1) << "Available networks could not be determined" << std::endl;

    struct ifaddrs *address = addresses;
    while(address) {
        if (address->ifa_addr != NULL && address->ifa_addr->sa_family == AF_INET) {
            char ap[100];
            const int family_size = sizeof(struct sockaddr_in);
            getnameinfo(address->ifa_addr,family_size, ap, sizeof(ap), 0, 0, NI_NUMERICHOST);
            networks[address->ifa_name] = ap;
        }
        address = address->ifa_next;
    }
    freeifaddrs(addresses);
#else
    CO_SIM_IO_ERROR << "GetIpv4Addresses is not supported for this operating system!" << std::endl;
#endif

    return networks;

    CO_SIM_IO_CATCH
}

std::string GetIpAddress(const Info& I_Settings)
{
    CO_SIM_IO_TRY

    if (I_Settings.Has("ip_address")) {
        return I_Settings.Get<std::string>("ip_address");
    } else if (I_Settings.Has("network_name")) {
        const std::string requested_network = I_Settings.Get<std::string>("network_name");
        const auto avail_networks = GetIpv4Addresses();
        const auto it_ip_address = avail_networks.find(requested_network);

        if (it_ip_address != avail_networks.end()) {
            return it_ip_address->second;
        } else {
            std::stringstream err_msg;
            err_msg << "The network with name \"" << requested_network << "\" could not be found! Only the following networks are available:";
            for (const auto& r_entr : avail_networks) {
                err_msg << "\n    Network name: " << r_entr.first << " | IP address: " << r_entr.second;
            }
            CO_SIM_IO_ERROR << err_msg.str() << std::endl;
        }
    } else {
        return LOCAL_IP_ADDRESS; // local loopback interface
    }

    CO_SIM_IO_CATCH
}

struct ConnectionInfo
{
    int PortNumber;
    std::string IpAddress;

    friend class CoSimIO::Internals::Serializer;

    void save(CoSimIO::Internals::Serializer& rSerializer) const
    {
        rSerializer.save("PortNumber", PortNumber);
        rSerializer.save("IpAddress", IpAddress);
    }

    void load(CoSimIO::Internals::Serializer& rSerializer)
    {
        rSerializer.load("PortNumber", PortNumber);
        rSerializer.load("IpAddress", IpAddress);
    }
};

// this is required as the serializer cannot handle newlines
void PrepareStringForAsciiSerialization(std::string& rString)
{
    const char disallowed_chars[] = {'<', '>'};
    for (const auto ch : disallowed_chars) {
        CO_SIM_IO_ERROR_IF_NOT(rString.find(ch) == std::string::npos) << "String contains a character that is not allowed: \"" << std::string(1,ch) << "\"!" << std::endl;
    }

    std::replace(rString.begin(), rString.end(), '\n', '<');
    std::replace(rString.begin(), rString.end(), '"', '>');
}

void RevertAsciiSerialization(std::string& rString)
{
    std::replace(rString.begin(), rString.end(), '<', '\n');
    std::replace(rString.begin(), rString.end(), '>', '"');
}

} // helpers namespace

SocketCommunication::SocketCommunication(
    const Info& I_Settings,
    std::shared_ptr<DataCommunicator> I_DataComm)
    : BaseType(I_Settings, I_DataComm)
{
    CO_SIM_IO_TRY

    if (GetIsPrimaryConnection()) {
        mIpAddress = GetIpAddress(I_Settings);
    }

    CO_SIM_IO_CATCH
}

Info SocketCommunication::ConnectDetail(const Info& I_Info)
{
    CO_SIM_IO_TRY

    if (!GetIsPrimaryConnection()) {GetConnectionInformation();}

    CO_SIM_IO_INFO_IF("CoSimIO", GetDataCommunicator().IsDistributed() && GetDataCommunicator().Rank()==0 && mIpAddress==LOCAL_IP_ADDRESS) << "Warning: Using the local IP address when connecting with MPI, this does not work in a distributed memory machine when communicating between different compute nodes!\nEither directly specify the IP address (with \"ip_address\") or specify the name of the network to be used (with \"network_name\")!" << std::endl;

    CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Using IP-Address: " << mIpAddress << " and port number: " << mPortNumber << std::endl;

    using namespace asio::ip;

    mpAsioSocket = std::make_shared<asio::ip::tcp::socket>(mAsioContext);
    if (GetIsPrimaryConnection()) { // this is the server
        mpAsioAcceptor->accept(*mpAsioSocket);
        mpAsioAcceptor->close();
        mpAsioAcceptor.reset();
    } else { // this is the client
        tcp::endpoint my_endpoint(asio::ip::make_address(mIpAddress), mPortNumber);
        mpAsioSocket->connect(my_endpoint);
    }

    return BaseType::ConnectDetail(I_Info);

    CO_SIM_IO_CATCH
}

void SocketCommunication::PrepareConnection(const Info& I_Info)
{
    CO_SIM_IO_TRY

    // preparing the acceptors to get the ports used for connecting the sockets
    if (GetIsPrimaryConnection()) {
        using namespace asio::ip;
        tcp::endpoint port_selection_endpoint(asio::ip::make_address(mIpAddress), 0); // using port 0 means that it will look for a free port
        mpAsioAcceptor = std::make_shared<tcp::acceptor>(mAsioContext, port_selection_endpoint);
        mPortNumber = mpAsioAcceptor->local_endpoint().port();

        // collect all IP-addresses and port numbers on rank 0 to
        // exchange them during the handshake (which happens only on rank 0)

        const auto& r_data_comm = GetDataCommunicator();
        std::vector<ConnectionInfo> conn_infos;

        ConnectionInfo my_conn_info {mPortNumber, mIpAddress};

        if (r_data_comm.Rank() == 0) {
            conn_infos.resize(r_data_comm.Size());
            conn_infos[0] = {mPortNumber, mIpAddress};
            for (int i=1; i<r_data_comm.Size(); ++i) {
                r_data_comm.Recv(conn_infos[i], i);
            }
        } else {
            r_data_comm.Send(my_conn_info, 0);
        }

        StreamSerializer serializer(GetSerializerTraceType());
        serializer.save("conn_info", conn_infos);
        mSerializedConnectionInfo = serializer.GetStringRepresentation();
        if (GetSerializerTraceType() != Serializer::TraceType::SERIALIZER_NO_TRACE) {
            PrepareStringForAsciiSerialization(mSerializedConnectionInfo);
        }
    }

    CO_SIM_IO_CATCH
}

Info SocketCommunication::GetCommunicationSettings() const
{
    CO_SIM_IO_TRY

    Info info;

    if (GetIsPrimaryConnection() && GetDataCommunicator().Rank() == 0) {
        info.Set("connection_info", mSerializedConnectionInfo);
    }

    return info;

    CO_SIM_IO_CATCH
}

void SocketCommunication::GetConnectionInformation()
{
    CO_SIM_IO_TRY

    CO_SIM_IO_ERROR_IF(GetIsPrimaryConnection()) << "This function can only be used as secondary connection!" << std::endl;

    const auto partner_info = GetPartnerInfo();
    std::string serialized_info = partner_info.Get<Info>("communication_settings").Get<std::string>("connection_info");

    if (GetSerializerTraceType() != Serializer::TraceType::SERIALIZER_NO_TRACE) {
        RevertAsciiSerialization(serialized_info);
    }

    std::vector<ConnectionInfo> conn_infos;

    StreamSerializer serializer(serialized_info, GetSerializerTraceType());
    serializer.load("conn_info", conn_infos);

    CO_SIM_IO_ERROR_IF(static_cast<int>(conn_infos.size()) != GetDataCommunicator().Size()) << "Wrong number of connection infos!" << std::endl;

    const auto& my_conn_info = conn_infos[GetDataCommunicator().Rank()];
    mPortNumber = my_conn_info.PortNumber;
    mIpAddress = my_conn_info.IpAddress;

    CO_SIM_IO_CATCH
}

} // namespace Internals
} // namespace CoSimIO
