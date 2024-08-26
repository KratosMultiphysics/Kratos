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

std::unique_ptr<IntegrationScheme> MakeLobattoIntegrationScheme(std::size_t NumberOfPoints)
{
    return std::make_unique<LobattoIntegrationScheme>(NumberOfPoints);
}

double SumOfWeights(const Geo::IntegrationPointVectorType& rIntegrationPoints)
{
    auto weights = std::vector<double>{};
    std::transform(rIntegrationPoints.begin(), rIntegrationPoints.end(), std::back_inserter(weights),
                   [](const auto& rIntegrationPoint) { return rIntegrationPoint.Weight(); });
    return std::accumulate(weights.cbegin(), weights.cend(), 0.0);
}

void ExpectLocalCoordinatesIncludeRangeBounds(const Geo::IntegrationPointVectorType& rIntegrationPoints)
{
    constexpr auto tolerance = 1.0e-6;

    auto contains_lower_bound_of_xi = [tolerance](const auto& rPoint) {
        return std::abs(rPoint[0] + 1.0) <= tolerance;
    };
    KRATOS_EXPECT_TRUE(std::any_of(rIntegrationPoints.begin(), rIntegrationPoints.end(), contains_lower_bound_of_xi))

    auto contains_upper_bound_of_xi = [tolerance](const auto& rPoint) {
        return std::abs(rPoint[0] - 1.0) <= tolerance;
    };
    KRATOS_EXPECT_TRUE(std::any_of(rIntegrationPoints.begin(), rIntegrationPoints.end(), contains_upper_bound_of_xi))
}

void ExpectLocalCoordinatesAreInRange(const Geo::IntegrationPointVectorType& rIntegrationPoints)
{
    constexpr auto tolerance = 1.0e-6;

    auto xi_is_in_range = [tolerance](const auto& rPoint) {
        return std::abs(rPoint[0]) - 1.0 <= tolerance;
    };
    KRATOS_EXPECT_TRUE(std::all_of(rIntegrationPoints.begin(), rIntegrationPoints.end(), xi_is_in_range))

    auto non_xi_coordinates_must_be_near_zero = [tolerance](const auto& rPoint) {
        return (std::abs(rPoint[1]) <= tolerance) && (std::abs(rPoint[2]) <= tolerance);
    };
    KRATOS_EXPECT_TRUE(std::all_of(rIntegrationPoints.begin(), rIntegrationPoints.end(),
                                   non_xi_coordinates_must_be_near_zero))
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(ALobattoIntegrationSchemeIsAnIntegrationScheme, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    constexpr auto number_of_points = 2;
    const auto     scheme           = LobattoIntegrationScheme{number_of_points};

    KRATOS_EXPECT_NE(dynamic_cast<const IntegrationScheme*>(&scheme), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(NumberOfIntegrationPointsMatchesTheNumberOfPointsGivenAtConstructionTime,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto supported_numbers_of_points = std::vector<std::size_t>{2, 3};

    for (auto number : supported_numbers_of_points) {
        const auto scheme = MakeLobattoIntegrationScheme(number);

        KRATOS_EXPECT_EQ(scheme->GetNumberOfIntegrationPoints(), number);
        KRATOS_EXPECT_EQ(scheme->GetIntegrationPoints().size(), number);
    }
}

KRATOS_TEST_CASE_IN_SUITE(SumOfIntegrationPointWeightsOfAllSupportedLobattoSchemesEqualsTwo,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto supported_numbers_of_points = std::vector<std::size_t>{2, 3};

    constexpr auto relative_tolerance = 1.0e-6;
    for (auto number : supported_numbers_of_points) {
        const auto scheme = MakeLobattoIntegrationScheme(number);

        KRATOS_EXPECT_RELATIVE_NEAR(SumOfWeights(scheme->GetIntegrationPoints()), 2.0, relative_tolerance)
    }
}

KRATOS_TEST_CASE_IN_SUITE(PointsOfAllSupportedLobattoSchemesMustBeInRangeAndIncludeBounds,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto supported_numbers_of_points = std::vector<std::size_t>{2, 3};

    for (auto number : supported_numbers_of_points) {
        const auto scheme = MakeLobattoIntegrationScheme(number);

        ExpectLocalCoordinatesAreInRange(scheme->GetIntegrationPoints());
        ExpectLocalCoordinatesIncludeRangeBounds(scheme->GetIntegrationPoints());
    }
}

KRATOS_TEST_CASE_IN_SUITE(AttemptingToMakeALobattoSchemeWithAnUnsupportedNumberOfPointsRaisesAnError,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto some_unsupported_numbers_of_points = std::vector<std::size_t>{0, 1, 4, 5, 6, 7};

    for (auto number : some_unsupported_numbers_of_points) {
        const auto expected_error_message =
            "Can't construct Lobatto integration scheme: no support for " + std::to_string(number) + " point(s)";
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(LobattoIntegrationScheme{number}, expected_error_message)
    }
}

} // namespace Kratos::Testing