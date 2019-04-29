#pragma once

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

constexpr int GetMajorVersion() {
    return KRATOS_MAJOR_VERSION;
}

constexpr int GetMinorVersion() {
    return KRATOS_MINOR_VERSION;
}

constexpr int GetPatchVersion() {
    return KRATOS_PATCH_VERSION;
}

const char * GetCommit();
const char * GetBuildType();
const char * GetVersionString();

} // namespace Kratos