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

namespace
{

void ExpectIntegrationPointsAreNear(const Geo::IntegrationPointVectorType& rExpectedIntegrationPoints,
                                    const Geo::IntegrationPointVectorType& rActualIntegrationPoints,
                                    double                                 RelativeTolerance)
{
    KRATOS_EXPECT_EQ(rExpectedIntegrationPoints.size(), rActualIntegrationPoints.size());

    for (auto i = 0; i < rExpectedIntegrationPoints.size(); ++i) {
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(rExpectedIntegrationPoints[i], rActualIntegrationPoints[i], RelativeTolerance)
        KRATOS_EXPECT_RELATIVE_NEAR(rExpectedIntegrationPoints[i].Weight(),
                                    rActualIntegrationPoints[i].Weight(), RelativeTolerance)
    }
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(ALobattoIntegrationSchemeIsAnIntegrationScheme, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto lobatto_integration_scheme = LobattoIntegrationScheme{};

    KRATOS_EXPECT_NE(dynamic_cast<const IntegrationScheme*>(&lobatto_integration_scheme), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(ADefaultConstructedLobattoIntegrationSchemeHasNoIntegrationPoints,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto lobatto_integration_scheme = LobattoIntegrationScheme{};

    KRATOS_EXPECT_EQ(lobatto_integration_scheme.GetNumberOfIntegrationPoints(), 0);
    KRATOS_EXPECT_TRUE(lobatto_integration_scheme.GetIntegrationPoints().empty())
}

KRATOS_TEST_CASE_IN_SUITE(CantConstructALobattoIntegrationSchemeWithLessThanTwoPoints, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(LobattoIntegrationScheme{0},
                                      "Can't construct Lobatto integration scheme: got 0 point(s)")
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(LobattoIntegrationScheme{1},
                                      "Can't construct Lobatto integration scheme: got 1 point(s)")
}

KRATOS_TEST_CASE_IN_SUITE(ALobattoIntegrationSchemeConstructedFromIntegrationPointsCanReturnThem,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto integration_points = Geo::IntegrationPointVectorType{{-1.0, 1.0}, {1.0, 1.0}};
    const auto lobatto_integration_scheme = LobattoIntegrationScheme{integration_points};

    KRATOS_EXPECT_EQ(lobatto_integration_scheme.GetNumberOfIntegrationPoints(), integration_points.size());

    constexpr auto relative_tolerance = 1.0e-6;
    ExpectIntegrationPointsAreNear(
        integration_points, lobatto_integration_scheme.GetIntegrationPoints(), relative_tolerance);
}

} // namespace Kratos::Testing