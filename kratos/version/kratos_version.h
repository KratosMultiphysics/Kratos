#pragma once

#include <string>

class KratosVersion {
    KratosVersion() = delete;
    static constexpr int GetMajor();
    static constexpr int GetMinor();
    static constexpr int GetPatch();
    static constexpr const char * GetCommit();
    static constexpr const char * GetBuildType();
    static constexpr const char * GetVersionString();
};