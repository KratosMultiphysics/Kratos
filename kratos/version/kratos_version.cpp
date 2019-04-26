#include "version/kratos_version.h"

// Kratos Minor, major and patch
#ifndef KRATOS_MAJOR_VERSION
#define KRATOS_MAJOR_VERSION 7
#endif

#ifndef KRATOS_MINOR_VERSION
#define KRATOS_MINOR_VERSION 0
#endif

#ifndef KRATOS_PATCH_VERSION
#define KRATOS_PATCH_VERSION 0
#endif

// GiT revision at configure
#ifndef KRATOS_SHA1_NUMBER
#define KRATOS_SHA1_NUMBER "0"
#endif

// Build type
#ifndef KRATOS_BUILD_TYPE
#define KRATOS_BUILD_TYPE "Release"
#endif

// Full version
#define KRATOS_TO_STRING(X) #X
#define KRATOS_VERSION \
KRATOS_TO_STRING(KRATOS_MAJOR_VERSION) "." \
KRATOS_TO_STRING(KRATOS_MINOR_VERSION) "." \
KRATOS_TO_STRING(KRATOS_PATCH_VERSION) "-" \
KRATOS_SHA1_NUMBER "-" \
KRATOS_BUILD_TYPE

constexpr int KratosVersion::GetMajor() {
    return KRATOS_MAJOR_VERSION;
}

constexpr int GetMinor() {
    return KRATOS_MINOR_VERSION;
}

constexpr int GetPatch() {
    return KRATOS_PATCH_VERSION;
}

constexpr const char * GetCommit() {
    return KRATOS_SHA1_NUMBER;
}

constexpr const char * GetBuildType() {
    return KRATOS_BUILD_TYPE;
}

constexpr const char * GetVersionString() {
    return KRATOS_VERSION;
}