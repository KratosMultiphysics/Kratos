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

// Project includes
#include "includes/kratos_version.h"

namespace Kratos {

// Kratos patch
#ifndef KRATOS_PATCH_VERSION
#define KRATOS_PATCH_VERSION "0"
#endif

// GiT revision at configure
#ifndef KRATOS_SHA1_NUMBER
#define KRATOS_SHA1_NUMBER "0"
#endif

// Build type
#ifndef KRATOS_BUILD_TYPE
#define KRATOS_BUILD_TYPE "Release"
#endif

// Architecture type
#if defined(__x86_64__) || defined(_M_X64) || defined(__amd64)
#define KRATOS_ARCH_TYPE "x86_64"
#elif defined(__i386__) || defined(_M_IX86) || defined(_X86_) || defined(__i386)
#define KRATOS_ARCH_TYPE "x86"
#elif defined(__arm__)
#define KRATOS_ARCH_TYPE "ARM32"
#elif defined(__aarch64__)
#define KRATOS_ARCH_TYPE "ARM64"
#else
#define KRATOS_ARCH_TYPE "Unknown architecture"
#endif

// Full version
#define KRATOS_TO_STRING_(X) #X
#define KRATOS_TO_STRING(X) KRATOS_TO_STRING_(X)
#define KRATOS_VERSION_STRING \
KRATOS_TO_STRING(KRATOS_MAJOR_VERSION) "." \
KRATOS_TO_STRING(KRATOS_MINOR_VERSION) "." \
KRATOS_TO_STRING(KRATOS_PATCH_VERSION) "-" \
KRATOS_SHA1_NUMBER "-" \
KRATOS_BUILD_TYPE  "-" \
KRATOS_ARCH_TYPE

std::string GetPatchVersion() {
    return KRATOS_PATCH_VERSION;
}

std::string GetCommitVersion() {
    return KRATOS_SHA1_NUMBER;
}

std::string GetBuildType() {
    return KRATOS_BUILD_TYPE;
}

std::string GetVersionString() {
    return KRATOS_VERSION_STRING;
}

} // namespace Kratos
