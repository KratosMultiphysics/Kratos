#include "version/kratos_version.h"

namespace Kratos {

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

int GetMajorVersion() {
    return KRATOS_MAJOR_VERSION;
}

int GetMinorVersion() {
    return KRATOS_MINOR_VERSION;
}

int GetPatchVersion() {
    return KRATOS_PATCH_VERSION;
}

const char * GetCommitVersion() {
    return KRATOS_SHA1_NUMBER;
}

const char * GetBuildType() {
    return KRATOS_BUILD_TYPE;
}

const char * GetVersionString() {
    return KRATOS_VERSION;
}

} // namespace Kratos