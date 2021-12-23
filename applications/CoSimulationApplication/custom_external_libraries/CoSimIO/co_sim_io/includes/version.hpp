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

#ifndef CO_SIM_IO_VERSION_INCLUDED
#define CO_SIM_IO_VERSION_INCLUDED

// System includes
#include <string>

// Project includes
#include "co_sim_io_api.hpp"

namespace CoSimIO {

int CO_SIM_IO_API GetMajorVersion();

int CO_SIM_IO_API GetMinorVersion();

std::string CO_SIM_IO_API GetPatchVersion();

} // namespace CoSimIO

#endif // CO_SIM_IO_VERSION_INCLUDED
