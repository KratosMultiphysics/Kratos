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
#include "custom_elements/lumped_integration_scheme.h"
#include "includes/expect.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

#include <algorithm>
#include <cstdlib>
#include <memory>
#include <numeric>
#include <string>

using namespace Kratos;

namespace
{

std::unique_ptr<IntegrationScheme> MakeLobattoIntegrationScheme(std::size_t NumberOfPoints)
{
    return std::make_unique<LobattoIntegrationScheme>(NumberOfPoints);
}

std::vector<std::size_t> SupportedNumbersOfPointsForLobattoIntegration() { return {2, 3}; }

std::unique_ptr<IntegrationScheme> MakeLumpedIntegrationScheme(std::size_t NumberOfPoints)
{
    return std::make_unique<LumpedIntegrationScheme>(NumberOfPoints);
}

std::vector<std::size_t> SupportedNumbersOfPointsForTriangleLumpedIntegration() { return {3, 6}; }
std::vector<std::size_t> SupportedNumbersOfPointsForQuadrilateralLumpedIntegration() { return {4, 8}; }

std::vector<std::size_t> SupportedNumbersOfPointsForLumpedIntegration() {
    auto sup_triang_num = SupportedNumbersOfPointsForTriangleLumpedIntegration();
    auto sup_quad_num   = SupportedNumbersOfPointsForQuadrilateralLumpedIntegration();
    auto sup_num = sup_triang_num;
    std::move(sup_quad_num.begin(), sup_quad_num.end(), std::back_inserter(sup_num));
    return sup_num;
}

double SumOfWeights(const Geo::IntegrationPointVectorType& rIntegrationPoints)
{
    auto weights = std::vector<double>{};
    weights.reserve(rIntegrationPoints.size());
    std::transform(rIntegrationPoints.begin(), rIntegrationPoints.end(), std::back_inserter(weights),
                   [](const auto& rIntegrationPoint) { return rIntegrationPoint.Weight(); });
    return std::accumulate(weights.cbegin(), weights.cend(), 0.0);
}

void ExpectLocalCoordinatesIncludeRangeBounds(const Geo::IntegrationPointVectorType& rIntegrationPoints, const size_t direction_index)
{
    constexpr auto tolerance = 1.0e-6;

    auto is_at_lower_bound_of_xi = [tolerance, direction_index](const auto& rPoint) {
        return std::abs(rPoint[direction_index] + 1.0) <= tolerance;
    };
    KRATOS_EXPECT_TRUE(std::any_of(rIntegrationPoints.begin(), rIntegrationPoints.end(), is_at_lower_bound_of_xi))

    auto is_at_upper_bound_of_xi = [tolerance, direction_index](const auto& rPoint) {
        return std::abs(rPoint[direction_index] - 1.0) <= tolerance;
    };
    KRATOS_EXPECT_TRUE(std::any_of(rIntegrationPoints.begin(), rIntegrationPoints.end(), is_at_upper_bound_of_xi))
}

void ExpectLineLocalCoordinatesAreInRange(const Geo::IntegrationPointVectorType& rIntegrationPoints)
{
    constexpr auto tolerance = 1.0e-6;

    auto xi_is_in_range = [tolerance](const auto& rPoint) {
        return std::abs(rPoint[0]) - 1.0 <= tolerance;
    };
    KRATOS_EXPECT_TRUE(std::all_of(rIntegrationPoints.begin(), rIntegrationPoints.end(), xi_is_in_range))
    // The isoparametric space is 1D (xi only), so the other parameters should be equal to zero
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
    for (auto number : SupportedNumbersOfPointsForLobattoIntegration()) {
        const auto scheme = MakeLobattoIntegrationScheme(number);

        KRATOS_EXPECT_EQ(scheme->GetNumberOfIntegrationPoints(), number);
        KRATOS_EXPECT_EQ(scheme->GetIntegrationPoints().size(), number);
    }
}

KRATOS_TEST_CASE_IN_SUITE(SumOfIntegrationPointWeightsOfAllSupportedLobattoSchemesEqualsTwo,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    constexpr auto relative_tolerance = 1.0e-6;
    for (auto number : SupportedNumbersOfPointsForLobattoIntegration()) {
        const auto scheme = MakeLobattoIntegrationScheme(number);

        KRATOS_EXPECT_RELATIVE_NEAR(SumOfWeights(scheme->GetIntegrationPoints()), 2.0, relative_tolerance)
    }
}

KRATOS_TEST_CASE_IN_SUITE(PointsOfAllSupportedLobattoSchemesMustBeInRangeAndIncludeBounds,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    for (auto number : SupportedNumbersOfPointsForLobattoIntegration()) {
        const auto scheme = MakeLobattoIntegrationScheme(number);

        ExpectLineLocalCoordinatesAreInRange(scheme->GetIntegrationPoints());
        ExpectLocalCoordinatesIncludeRangeBounds(scheme->GetIntegrationPoints(), size_t(0));
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
KRATOS_TEST_CASE_IN_SUITE(ALumpedIntegrationSchemeIsAnIntegrationScheme, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    constexpr auto number_of_points = 3;
    const auto     scheme           = LumpedIntegrationScheme{number_of_points};

    KRATOS_EXPECT_NE(dynamic_cast<const IntegrationScheme*>(&scheme), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(NumberOfLumpedIntegrationPointsMatchesTheNumberOfPointsGivenAtConstructionTime,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    for (auto number : SupportedNumbersOfPointsForLumpedIntegration()) {
        const auto scheme = MakeLumpedIntegrationScheme(number);

        KRATOS_EXPECT_EQ(scheme->GetNumberOfIntegrationPoints(), number);
        KRATOS_EXPECT_EQ(scheme->GetIntegrationPoints().size(), number);
    }
}

KRATOS_TEST_CASE_IN_SUITE(SumOfIntegrationPointWeightsOfAllSupportedLumpedSchemesEqualsOne,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    constexpr auto relative_tolerance = 1.0e-6;
    for (auto number : SupportedNumbersOfPointsForLumpedIntegration()) {
        const auto scheme = MakeLumpedIntegrationScheme(number);

        KRATOS_EXPECT_RELATIVE_NEAR(SumOfWeights(scheme->GetIntegrationPoints()), 1.0, relative_tolerance)
    }
}

KRATOS_TEST_CASE_IN_SUITE(AttemptingToMakeALumpedSchemeWithAnUnsupportedNumberOfPointsRaisesAnError,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto some_unsupported_numbers_of_points = std::vector<std::size_t>{0, 1, 2, 5, 7, 9};

    for (auto number : some_unsupported_numbers_of_points) {
        const auto expected_error_message =
            "Can't construct Lumped integration scheme: no support for " + std::to_string(number) + " point(s)";
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(LumpedIntegrationScheme{number}, expected_error_message)
    }
}

KRATOS_TEST_CASE_IN_SUITE(PointsOfAllSupportedQuadrilateralLumpedSchemesMustBeInRangeAndIncludeBounds,
                              KratosGeoMechanicsFastSuiteWithoutKernel)
{
    for (auto number : SupportedNumbersOfPointsForQuadrilateralLumpedIntegration()) {
        const auto scheme = MakeLumpedIntegrationScheme(number);

        ExpectLocalCoordinatesIncludeRangeBounds(scheme->GetIntegrationPoints(), size_t(0));
        ExpectLocalCoordinatesIncludeRangeBounds(scheme->GetIntegrationPoints(), size_t(1));
    }
}
} // namespace Kratos::Testing