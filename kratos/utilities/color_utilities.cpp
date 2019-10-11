//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
// #ifdef KRATOS_COMPILED_IN_WINDOWS
    #include <stdio.h>
    #ifndef _WIN32_WINNT
        #define _WIN32_WINNT 0x0601 // Windows 7 or later. See https://docs.microsoft.com/en-us/cpp/porting/modifying-winver-and-win32-winnt?view=vs-2019
    #endif
    #include <windows.h>
    #if defined(__MINGW32__) || defined(__MINGW64__)
        #include <winbase.h>
    #endif
// #endif

// External includes

// Project includes
#include "utilities/color_utilities.h"

namespace Kratos
{
namespace ColorUtilities
{
void InitializeSTDOutput()
{
#ifdef KRATOS_COMPILED_IN_WINDOWS
    SetConsoleMode(GetStdHandle(STD_OUTPUT_HANDLE), ENABLE_VIRTUAL_TERMINAL_PROCESSING);
#endif
}

} // namespace ColorUtilities
} // namespace Kratos
