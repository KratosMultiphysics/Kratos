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

double SumOfWeights(const Geo::IntegrationPointVectorType& rIntegrationPoints)
{
    auto weights = std::vector<double>{};
    std::transform(rIntegrationPoints.begin(), rIntegrationPoints.end(), std::back_inserter(weights),
                   [](const auto& rIntegrationPoint) { return rIntegrationPoint.Weight(); });
    return std::accumulate(weights.cbegin(), weights.cend(), 0.0);
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(ALobattoIntegrationSchemeIsAnIntegrationScheme, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto lobatto_integration_scheme = LobattoIntegrationScheme{};

    KRATOS_EXPECT_NE(dynamic_cast<const IntegrationScheme*>(&lobatto_integration_scheme), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(ADefaultConstructedLobattoIntegrationSchemeHasTwoIntegrationPoints,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto lobatto_integration_scheme = LobattoIntegrationScheme{};

    KRATOS_EXPECT_EQ(lobatto_integration_scheme.GetNumberOfIntegrationPoints(), 2);
}

KRATOS_TEST_CASE_IN_SUITE(CantConstructALobattoIntegrationSchemeWhenNumberOfPointsIsNotEqualToTwo,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto some_unsupported_numbers_of_nodes = std::vector<std::size_t>{0, 1, 3, 4, 5, 6, 7};

    for (auto number : some_unsupported_numbers_of_nodes) {
        const auto expected_error_message =
            "Can't construct Lobatto integration scheme: no support for " + std::to_string(number) + " point(s)";
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(LobattoIntegrationScheme{number}, expected_error_message)
    }
}

KRATOS_TEST_CASE_IN_SUITE(ATwoPointLobattoIntegrationSchemeUsesEndPointsOfLineWithUnityWeight,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto lobatto_integration_scheme = LobattoIntegrationScheme{2};

    const auto expected_integration_points = Geo::IntegrationPointVectorType{{-1.0, 1.0}, {1.0, 1.0}};
    constexpr auto relative_tolerance = 1.0e-6;
    ExpectIntegrationPointsAreNear(expected_integration_points,
                                   lobatto_integration_scheme.GetIntegrationPoints(), relative_tolerance);

    KRATOS_EXPECT_RELATIVE_NEAR(SumOfWeights(lobatto_integration_scheme.GetIntegrationPoints()), 2.0, relative_tolerance)
}

} // namespace Kratos::Testing