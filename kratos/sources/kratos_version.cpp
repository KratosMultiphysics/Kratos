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

// GiT branch name at configure
#ifndef KRATOS_BRANCH_NAME
#define KRATOS_BRANCH_NAME ""
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
#ifndef KRATOS_TO_STRING_
    #define KRATOS_TO_STRING_(X) #X
#endif
#ifndef KRATOS_TO_STRING
    #define KRATOS_TO_STRING(X) KRATOS_TO_STRING_(X)
#endif
#define KRATOS_VERSION_STRING \
KRATOS_TO_STRING(KRATOS_MAJOR_VERSION) "." \
KRATOS_TO_STRING(KRATOS_MINOR_VERSION) "." \
KRATOS_TO_STRING(KRATOS_PATCH_VERSION) "-" \
KRATOS_BRANCH_NAME "-" \
KRATOS_SHA1_NUMBER "-" \
KRATOS_BUILD_TYPE  "-" \
KRATOS_ARCH_TYPE

// Define OS name
#if defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
    #define KRATOS_OS_NAME "GNU/Linux"
#elif defined(__APPLE__) && defined(__MACH__)
    #define KRATOS_OS_NAME "Mac OS"
#elif defined(_WIN32) || defined(_WIN64)
    #define KRATOS_OS_NAME "Windows"
#else
    #define KRATOS_OS_NAME "Unknown OS"
#endif

// Define compiler label
#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
    // the intel compiler internally uses gcc or clang
    // therefore must check it first, otherwise it detects wrong compiler
    #define KRATOS_COMPILER_LABEL "Intel-" \
    KRATOS_TO_STRING(__VERSION__)
#elif defined(__clang__)
    #define KRATOS_COMPILER_LABEL "Clang-" \
    KRATOS_TO_STRING(__clang_major__) \
    "." \
    KRATOS_TO_STRING(__clang_minor__)
#elif defined(__GNUC__) || defined(__GNUG__)
    #define KRATOS_COMPILER_LABEL "GCC-" \
    KRATOS_TO_STRING(__GNUC__) \
    "." \
    KRATOS_TO_STRING(__GNUC_MINOR__)
#elif defined(_MSC_VER)
    #define KRATOS_COMPILER_LABEL "MSVC-" \
    KRATOS_TO_STRING(_MSC_VER)
#else
    #define KRATOS_COMPILER_LABEL "Unknown compiler"
#endif

std::string GetPatchVersion() {
    return KRATOS_PATCH_VERSION;
}

std::string GetCommit() {
    return KRATOS_SHA1_NUMBER;
}

std::string GetBranchName() {
    return KRATOS_BRANCH_NAME;
}

std::string GetBuildType() {
    return KRATOS_BUILD_TYPE;
}

std::string GetVersionString() {
    return KRATOS_VERSION_STRING;
}

std::string GetOSName() {
    return KRATOS_OS_NAME;
}

std::string GetCompiler() {
    return KRATOS_COMPILER_LABEL;
}

} // namespace Kratos
