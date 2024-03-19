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

#include <string>
#include <filesystem>

namespace Kratos::Testing
{

class TestUtilities
{
public:
    static bool CompareFiles(const std::filesystem::path& rPath1,
                             const std::filesystem::path& rPath2);
};

}
