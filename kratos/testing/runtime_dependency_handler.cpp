#include "runtime_dependency_handler.h"

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
        #ifdef _WIN32
        const std::string lib_file = rLibraryName+".dll";
        LibraryHandle library = LoadLibrary(TEXT(lib_file.c_str()));
        #else
        const std::string lib_file = "lib"+rLibraryName+".so";
        LibraryHandle library = dlopen(lib_file.c_str(), RTLD_NOW|RTLD_GLOBAL);
        #endif
        KRATOS_ERROR_IF(library == nullptr) << "Could not load runtime dependency library " << rLibraryName << std::endl;
        mLoadedLibraries.insert_or_assign(rApplicationName, library);
    }
}


RuntimeDependencyHandler::ApplicationMap RuntimeDependencyHandler::CreateApplications()
{
    using KratosApplication = Kratos::KratosApplication;
    typedef KratosApplication* (*CreateApplicationPointer)();

    ApplicationMap applications;

    for (auto& r_loaded: mLoadedLibraries) {
        LibraryHandle library = r_loaded.second;
        #ifdef _WIN32
        auto p_create_app = reinterpret_cast<CreateApplicationPointer>(GetProcAddress(library, "CreateApplication"));
        #else
        auto p_create_app = reinterpret_cast<CreateApplicationPointer>(dlsym(library, "CreateApplication"));
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
        #ifdef _WIN32
        FreeLibrary(iter->second);
        #else
        dlclose(iter->second);
        #endif
    }

    mLoadedLibraries.clear();
}

}