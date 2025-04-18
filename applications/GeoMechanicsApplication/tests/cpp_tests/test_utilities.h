// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#pragma once

#include <filesystem>

namespace Kratos::Testing
{

namespace Defaults
{

constexpr auto absolute_tolerance = 1.0e-12;
constexpr auto relative_tolerance = 1.0e-6;

} // namespace Defaults

class TestUtilities
{
public:
    static bool CompareFiles(const std::filesystem::path& rPath1, const std::filesystem::path& rPath2);
};

} // namespace Kratos::Testing
