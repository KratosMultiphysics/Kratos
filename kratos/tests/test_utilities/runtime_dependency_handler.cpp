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

#include "runtime_dependency_handler.h"

#include "includes/define.h" // defines KRATOS_COMPILED_IN... flags

#ifdef KRATOS_COMPILED_IN_WINDOWS
#include <windows.h>
#else
#include <dlfcn.h>
#endif

namespace Kratos::Testing {


RuntimeDependencyHandler::RuntimeDependencyHandler() {}


RuntimeDependencyHandler::~RuntimeDependencyHandler()
{
    ReleaseDependencies();
}


void RuntimeDependencyHandler::LoadDependency(const std::string& rApplicationName, const std::string& rLibraryName)
{
    auto found = mLoadedLibraries.find(rApplicationName);
    if (found == mLoadedLibraries.end()) {
        #ifdef KRATOS_COMPILED_IN_WINDOWS
        const std::string lib_file = rLibraryName+".dll";
        void* p_library = reinterpret_cast<void*>(LoadLibrary(TEXT(lib_file.c_str())));
        #else
        const std::string lib_file = "lib"+rLibraryName+".so";
        void* p_library = dlopen(lib_file.c_str(), RTLD_NOW|RTLD_GLOBAL);
        #endif
        KRATOS_ERROR_IF(p_library == nullptr) << "Could not load runtime dependency library " << rLibraryName << std::endl;
        mLoadedLibraries.insert({rApplicationName, p_library});
    }
}


RuntimeDependencyHandler::ApplicationMap RuntimeDependencyHandler::CreateApplications()
{
    using KratosApplication = Kratos::KratosApplication;
    typedef KratosApplication* (*CreateApplicationPointer)();

    ApplicationMap applications;

    for (auto& r_loaded: mLoadedLibraries) {
        void* p_library = r_loaded.second;
        #ifdef KRATOS_COMPILED_IN_WINDOWS
        auto p_create_app = reinterpret_cast<CreateApplicationPointer>(GetProcAddress(reinterpret_cast<HINSTANCE>(p_library), "CreateApplication"));
        #else
        auto p_create_app = reinterpret_cast<CreateApplicationPointer>(dlsym(p_library, "CreateApplication"));
        #endif
        KRATOS_ERROR_IF(p_create_app == nullptr) << "Could not find CreateApplication method for runtime dependency library " << r_loaded.first << std::endl;
        applications.insert({r_loaded.first, KratosApplication::Pointer((*p_create_app)())});
    }

    return applications;
}


void RuntimeDependencyHandler::ReleaseDependencies()
{
    for (auto iter = mLoadedLibraries.begin(); iter != mLoadedLibraries.end(); ++iter)
    {
        #ifdef KRATOS_COMPILED_IN_WINDOWS
        FreeLibrary(reinterpret_cast<HINSTANCE>(iter->second));
        #else
        dlclose(iter->second);
        #endif
    }

    mLoadedLibraries.clear();
}

}