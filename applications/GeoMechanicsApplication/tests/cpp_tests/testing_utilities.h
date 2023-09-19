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
#include <vector>

class TestingUtilities {
public:
    static void Expect_Equal(const std::vector<std::string>& rA, const std::vector<std::string>& rB);

};
