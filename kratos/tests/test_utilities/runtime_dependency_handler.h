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

#include "includes/kratos_application.h"

namespace Kratos::Testing {

class RuntimeDependencyHandler
{
public:

    using ApplicationMap = std::unordered_map<std::string, Kratos::KratosApplication::Pointer>;

    RuntimeDependencyHandler();

    ~RuntimeDependencyHandler();

    void LoadDependency(const std::string& rApplicationName, const std::string& rLibraryName);

    ApplicationMap CreateApplications();

    void ReleaseDependencies();

private:

    std::unordered_map<std::string, void*> mLoadedLibraries = {};

};

}