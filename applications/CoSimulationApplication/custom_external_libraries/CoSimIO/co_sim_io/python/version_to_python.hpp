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

#ifndef CO_SIM_IO_VERSION_TO_PYHON_INCLUDED
#define CO_SIM_IO_VERSION_TO_PYHON_INCLUDED

// Exposure of the CoSimIO to Python

// System includes
#include <string>
#include <sstream>

// pybind includes
#include <pybind11/pybind11.h>

// CoSimIO include
#include "co_sim_io.hpp"

namespace CoSimIO {

void AddCoSimIOVersionToPython(pybind11::module& m)
{
    namespace py = pybind11;

    std::stringstream version_stream;
    version_stream << CoSimIO::GetMajorVersion() << "."
                   << CoSimIO::GetMinorVersion() << "."
                   << CoSimIO::GetPatchVersion();
    m.attr("__version__") = version_stream.str();
}

} // namespace CoSimIO

#endif // CO_SIM_IO_VERSION_TO_PYHON_INCLUDED
