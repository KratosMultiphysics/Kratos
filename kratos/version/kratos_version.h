#pragma once

namespace Kratos {

int GetMajorVersion();
int GetMinorVersion();
int GetPatchVersion();
const char * GetCommit();
const char * GetBuildType();
const char * GetVersionString();

} // namespace Kratos