#include "version/kratos_version.h"

namespace Kratos {

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