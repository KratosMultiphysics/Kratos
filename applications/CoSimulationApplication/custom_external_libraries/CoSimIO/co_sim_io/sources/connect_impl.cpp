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
#include <string>
#include <unordered_map>

// Project includes
#include "includes/connect_impl.hpp"
#include "includes/utilities.hpp"

namespace CoSimIO {
namespace Internals {

static std::unordered_map<std::string, std::unique_ptr<Internals::Connection>> s_co_sim_connections;

bool HasConnection(const std::string& rConnectionName)
{
    return s_co_sim_connections.find(rConnectionName) != s_co_sim_connections.end();
}

Connection& GetConnection(const std::string& rConnectionName)
{
    CO_SIM_IO_ERROR_IF_NOT(HasConnection(rConnectionName)) << "Trying to use connection \"" << rConnectionName << "\" which does not exist!" << std::endl;
    return *s_co_sim_connections.at(rConnectionName);
}

void RemoveConnection(const std::string& rConnectionName)
{
    s_co_sim_connections.erase(rConnectionName);
}

Info ConnectImpl(const Info& I_Settings, std::shared_ptr<DataCommunicator> I_DataComm)
{
    const std::string my_name = I_Settings.Get<std::string>("my_name");
    const std::string connect_to = I_Settings.Get<std::string>("connect_to");

    // perform some checks
    Utilities::CheckEntry(my_name, "my_name");
    Utilities::CheckEntry(connect_to, "connect_to");
    CO_SIM_IO_ERROR_IF(my_name == connect_to) << "Connecting to self is not allowed!" << std::endl;

    const std::string connection_name = Utilities::CreateConnectionName(my_name, connect_to);

    CO_SIM_IO_ERROR_IF(HasConnection(connection_name)) << "A connection from \"" << my_name << "\" to \"" << connect_to << "\"already exists!" << std::endl;

    s_co_sim_connections[connection_name] = CoSimIO::make_unique<Connection>(I_Settings, I_DataComm);

    auto info = GetConnection(connection_name).Connect(I_Settings);
    info.Set<std::string>("connection_name", connection_name);

    return info;
}

} // namespace Internals
} // namespace CoSimIO
