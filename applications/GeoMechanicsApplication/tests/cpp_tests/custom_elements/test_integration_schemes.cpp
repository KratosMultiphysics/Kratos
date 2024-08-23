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

#include "custom_elements/lobatto_integration_scheme.h"
#include "includes/expect.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(ALobattoIntegrationSchemeIsAnIntegrationScheme, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto lobatto_integration_scheme = LobattoIntegrationScheme{};

    KRATOS_EXPECT_NE(dynamic_cast<const IntegrationScheme*>(&lobatto_integration_scheme), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(ADefaultConstructedLobattoIntegrationSchemeHasNoIntegrationPoints, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto lobatto_integration_scheme = LobattoIntegrationScheme{};

    KRATOS_EXPECT_EQ(lobatto_integration_scheme.GetNumberOfIntegrationPoints(), 0);
}

} // namespace Kratos::Testing