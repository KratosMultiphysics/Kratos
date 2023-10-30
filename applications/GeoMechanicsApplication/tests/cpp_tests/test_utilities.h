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

namespace Kratos::Testing
{

class TestUtilities
{
public:
    static bool CompareFiles(const std::string& p1, const std::string& p2);
};

}
