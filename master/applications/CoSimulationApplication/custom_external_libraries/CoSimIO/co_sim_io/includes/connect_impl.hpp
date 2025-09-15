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

#ifndef CO_SIM_IO_CONNECT_IMPL_INCLUDED
#define CO_SIM_IO_CONNECT_IMPL_INCLUDED

// System includes
#include <memory>
#include <iostream>

// Project includes
#include "define.hpp"
#include "info.hpp"
#include "connection.hpp"
#include "data_communicator.hpp"
#include "communication/factory.hpp"

namespace CoSimIO {
namespace Internals {

Info ConnectImpl(
    const Info& I_Settings,
    std::shared_ptr<DataCommunicator> I_DataComm,
    const CommunicationFactory& rCommFactory);

bool HasConnection(const std::string& rConnectionName);

Connection& GetConnection(const std::string& rConnectionName);

void RemoveConnection(const std::string& rConnectionName);

} // namespace Internals
} // namespace CoSimIO

#endif // CO_SIM_IO_CONNECT_IMPL_INCLUDED
