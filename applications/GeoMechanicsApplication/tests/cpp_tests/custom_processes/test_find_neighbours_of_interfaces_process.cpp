// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//

#include "custom_processes/find_neighbours_of_interfaces_process.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(FindNeighboursOfInterfacesProcess_IsAProcess, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto process = FindNeighboursOfInterfacesProcess{};
    EXPECT_NE(dynamic_cast<const Process*>(&process), nullptr);
}

} // namespace Kratos::Testing
