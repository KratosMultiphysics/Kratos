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

#ifndef CO_SIM_IO_VERSION_H_INCLUDED
#define CO_SIM_IO_VERSION_H_INCLUDED

// System includes
#include <string>

namespace CoSimIO {

constexpr int GetMajorVersion() {
    return 1;
}

constexpr int GetMinorVersion() {
    return 0;
}

std::string GetPatchVersion() {
    return "0";
}

} // namespace CoSimIO

#endif // CO_SIM_IO_VERSION_H_INCLUDED
