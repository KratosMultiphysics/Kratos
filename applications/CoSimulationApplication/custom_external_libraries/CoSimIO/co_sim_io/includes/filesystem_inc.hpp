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

#ifndef CO_SIM_IO_FILESYSTEM_INC_INCLUDED
#define CO_SIM_IO_FILESYSTEM_INC_INCLUDED

/* std::filesystem is part of C++17
While we use C++11, the alternative implementation from
"https://github.com/gulrak/filesystem" is used
*/

#ifndef WIN32_LEAN_AND_MEAN
    #define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
    #define NOMINMAX
#endif
#include "ghc/filesystem.hpp"
namespace fs = ghc::filesystem;

// use this once moving to C++17
// #include <filesystem>
// namespace fs = std::filesystem;

#endif // CO_SIM_IO_FILESYSTEM_INC_INCLUDED
