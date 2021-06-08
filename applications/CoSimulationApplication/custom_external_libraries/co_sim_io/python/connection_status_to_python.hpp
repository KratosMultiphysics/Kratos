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

#ifndef CO_SIM_IO_CONNECTION_STATUS_TO_PYHON_INCLUDED
#define CO_SIM_IO_CONNECTION_STATUS_TO_PYHON_INCLUDED

// Exposure of the CoSimIO to Python

// System includes

// pybind includes
#include <pybind11/pybind11.h>

// CoSimIO include
#include "../co_sim_io.hpp"

namespace CoSimIO {

void AddCoSimIOConnectionStatusToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::enum_<CoSimIO::ConnectionStatus>(m,"ConnectionStatus")
        .value("NotConnected",CoSimIO::ConnectionStatus::NotConnected)
        .value("Connected",CoSimIO::ConnectionStatus::Connected)
        .value("Disconnected", CoSimIO::ConnectionStatus::Disconnected)
        .value("ConnectionError", CoSimIO::ConnectionStatus::ConnectionError)
        .value("DisconnectionError", CoSimIO::ConnectionStatus::DisconnectionError)
        ;
}

} // namespace CoSimIO

#endif // CO_SIM_IO_CONNECTION_STATUS_TO_PYHON_INCLUDED
