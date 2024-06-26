//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela Dalmau
//
//

#pragma once

#include <unordered_map>

#include "includes/define.h"
#include "includes/kernel.h"

// Prevent Windows.h from polluting the global namespace
#ifdef KRATOS_COMPILED_IN_WINDOWS
    #define NOMINMAX
    #include <windows.h>
    #undef NOMINMAX
#else
    #include <dlfcn.h>
#endif

namespace Kratos::Testing {

class KRATOS_API(KRATOS_TEST_UTILS) RuntimeDependencyHandler
{
public:

    using ApplicationMap = std::unordered_map<std::string, Kratos::KratosApplication::Pointer>;

    RuntimeDependencyHandler();

    ~RuntimeDependencyHandler();

    void LoadDependency(const std::string& rApplicationName, const std::string& rLibraryName);

    ApplicationMap CreateApplications();

    void ReleaseDependencies();

private:

    #ifdef KRATOS_COMPILED_IN_WINDOWS
    using LibraryHandle = HINSTANCE;
    #else
    using LibraryHandle = void*;
    #endif

    std::unordered_map<std::string, LibraryHandle> mLoadedLibraries = {};
};

}