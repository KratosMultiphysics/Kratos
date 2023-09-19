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
#include "testing_utilities.h"
#include "testing/testing.h"

void TestingUtilities::Expect_Equal(const std::vector<std::string> &rA, const std::vector<std::string> &rB)
{
    KRATOS_EXPECT_EQ(rA.size(), rB.size());
    for (int i = 0; i < rA.size(); i++)
    {
        KRATOS_EXPECT_EQ(rA[i], rB[i]);
    }
}
