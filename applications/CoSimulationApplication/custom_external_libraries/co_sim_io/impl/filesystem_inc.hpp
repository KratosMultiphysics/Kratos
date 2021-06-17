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

/* This file selects which implementation of std::filesystem to be used
std::filesystem is part of C++17
When using only C++11 the alternative implementation from
"https://github.com/gulrak/filesystem" is used
*/

// To dynamically select std::filesystem where available, you could use:
#if defined(__cplusplus) && __cplusplus >= 201703L
    #if __has_include(<filesystem>) // has_include is C++17
        #define STD_FILESYSTEM_AVAILABLE
        #include <filesystem>
        namespace fs = std::filesystem;
    #endif
#endif
#ifndef STD_FILESYSTEM_AVAILABLE
    #define WIN32_LEAN_AND_MEAN
    #define NOMINMAX
    #include "../../external_libraries/ghc/filesystem.hpp"
    #undef NOMINMAX
    #undef WIN32_LEAN_AND_MEAN
    namespace fs = ghc::filesystem;
#endif

#undef STD_FILESYSTEM_AVAILABLE

#endif // CO_SIM_IO_FILESYSTEM_INC_INCLUDED
