//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/expect.h"
#include "utilities/geometry_tester.h"

namespace Kratos {
namespace Testing {

    KRATOS_TEST_CASE_IN_SUITE(GeometryTester, KratosCoreFastSuite)
    {
        Model this_model;
        KRATOS_EXPECT_TRUE(GeometryTesterUtility().RunTest(this_model));
    }

}  // namespace Testing.
}  // namespace Kratos.
