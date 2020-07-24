//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos Roig
//

// System includes
#include <string>

// External includes
#include "includes/kratos_export_api.h"

// Project includes

#pragma once

namespace Kratos {

// Kratos Minor and Major
#ifndef KRATOS_MAJOR_VERSION
#define KRATOS_MAJOR_VERSION 7
#endif

#ifndef KRATOS_MINOR_VERSION
#define KRATOS_MINOR_VERSION 1
#endif

constexpr int GetMajorVersion() {
    return KRATOS_MAJOR_VERSION;
}

constexpr int GetMinorVersion() {
    return KRATOS_MINOR_VERSION;
}

KRATOS_API(KRATOS_VERSION) std::string GetPatchVersion();
KRATOS_API(KRATOS_VERSION) std::string GetCommit();
KRATOS_API(KRATOS_VERSION) std::string GetBuildType();
KRATOS_API(KRATOS_VERSION) std::string GetVersionString();

} // namespace Kratos
