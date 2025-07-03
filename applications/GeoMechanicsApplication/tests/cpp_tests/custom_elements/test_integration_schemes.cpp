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

std::vector<std::size_t> SupportedNumbersOfPointsForQuadrilateralLumpedIntegration()
{
    return {4, 8};
}

std::vector<std::size_t> SupportedNumbersOfPointsForLumpedIntegration()
{
    auto sup_triang_num = SupportedNumbersOfPointsForTriangleLumpedIntegration();
    auto sup_quad_num   = SupportedNumbersOfPointsForQuadrilateralLumpedIntegration();
    auto sup_num        = sup_triang_num;
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

void ExpectLocalCoordinatesIncludeRangeBounds(const Geo::IntegrationPointVectorType& rIntegrationPoints,
                                              double      IsoLowerBound,
                                              double      IsoUpperBound,
                                              std::size_t DirectionIndex)
{
    constexpr auto tolerance = 1.0e-6;

    auto is_at_lower_bound_of_iso_coord = [tolerance, IsoLowerBound, DirectionIndex](const auto& rPoint) {
        // lower_bound is 0 or a negative number
        return std::abs(rPoint[DirectionIndex] - IsoLowerBound) <= tolerance;
    };
    KRATOS_EXPECT_TRUE(std::any_of(rIntegrationPoints.begin(), rIntegrationPoints.end(), is_at_lower_bound_of_iso_coord))

    auto is_at_upper_bound_of_iso_coord = [tolerance, IsoUpperBound, DirectionIndex](const auto& rPoint) {
        // upper_bound is 0 or a positive number
        return std::abs(rPoint[DirectionIndex] - IsoUpperBound) <= tolerance;
    };
    KRATOS_EXPECT_TRUE(std::any_of(rIntegrationPoints.begin(), rIntegrationPoints.end(), is_at_upper_bound_of_iso_coord))
}

void ExpectLocalCoordinatesAreInRange(const Geo::IntegrationPointVectorType& rIntegrationPoints,
                                      double                                 IsoLowerBound,
                                      double                                 IsoUpperBound,
                                      std::size_t                            DirectionIndex)
{
    constexpr auto tolerance = 1.0e-6;

    auto iso_coord_is_in_range = [tolerance, IsoLowerBound, IsoUpperBound, DirectionIndex](const auto& rPoint) {
        return rPoint[DirectionIndex] - IsoUpperBound <= tolerance &&
               IsoLowerBound - rPoint[DirectionIndex] <= tolerance;
    };
    KRATOS_EXPECT_TRUE(std::all_of(rIntegrationPoints.begin(), rIntegrationPoints.end(), iso_coord_is_in_range))
}

void ExpectLocalCoordinatesAreZero(const Geo::IntegrationPointVectorType& rIntegrationPoints, const size_t DirectionIndex)
{
    constexpr auto tolerance = 1.0e-6;

    auto iso_coordinate_must_be_near_zero = [tolerance, DirectionIndex](const auto& rPoint) {
        return (std::abs(rPoint[DirectionIndex]) <= tolerance);
    };
    KRATOS_EXPECT_TRUE(std::all_of(rIntegrationPoints.begin(), rIntegrationPoints.end(), iso_coordinate_must_be_near_zero))
}

std::vector<double> ComputeWeightsForLumpedIntegration(const size_t NumberOfPoints)
{
    const auto          scheme             = MakeLumpedIntegrationScheme(NumberOfPoints);
    auto                integration_points = scheme->GetIntegrationPoints();
    std::vector<double> weights;
    weights.reserve(integration_points.size());
    std::transform(integration_points.begin(), integration_points.end(), std::back_inserter(weights),
                   [](const auto& rIntegrationPoint) { return rIntegrationPoint.Weight(); });
    return weights;
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
    for (auto number : SupportedNumbersOfPointsForLobattoIntegration()) {
        const auto scheme = MakeLobattoIntegrationScheme(number);

        constexpr auto relative_tolerance = 1.0e-6;
        KRATOS_EXPECT_RELATIVE_NEAR(SumOfWeights(scheme->GetIntegrationPoints()), 2.0, relative_tolerance)
    }
}

KRATOS_TEST_CASE_IN_SUITE(PointsOfAllSupportedLobattoSchemesMustBeInRangeAndIncludeBounds,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    for (auto number : SupportedNumbersOfPointsForLobattoIntegration()) {
        const auto scheme = MakeLobattoIntegrationScheme(number);

        ExpectLocalCoordinatesAreInRange(scheme->GetIntegrationPoints(), -1.0, 1.0, std::size_t{0});
        ExpectLocalCoordinatesAreZero(scheme->GetIntegrationPoints(), std::size_t{1});
        ExpectLocalCoordinatesAreZero(scheme->GetIntegrationPoints(), std::size_t{2});
        ExpectLocalCoordinatesIncludeRangeBounds(scheme->GetIntegrationPoints(), -1.0, 1.0, std::size_t{0});
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
    for (auto number : SupportedNumbersOfPointsForLumpedIntegration()) {
        const auto scheme = MakeLumpedIntegrationScheme(number);

        constexpr auto relative_tolerance = 1.0e-6;
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

        ExpectLocalCoordinatesAreInRange(scheme->GetIntegrationPoints(), -1.0, 1.0, std::size_t{0});
        ExpectLocalCoordinatesAreInRange(scheme->GetIntegrationPoints(), -1.0, 1.0, std::size_t{1});
        ExpectLocalCoordinatesAreZero(scheme->GetIntegrationPoints(), std::size_t{2});
        ExpectLocalCoordinatesIncludeRangeBounds(scheme->GetIntegrationPoints(), -1.0, 1.0, std::size_t{0});
        ExpectLocalCoordinatesIncludeRangeBounds(scheme->GetIntegrationPoints(), -1.0, 1.0, std::size_t{1});
    }
}

KRATOS_TEST_CASE_IN_SUITE(PointsOfAllSupportedTriangleLumpedSchemesMustBeInRangeAndIncludeBounds,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    for (auto number : SupportedNumbersOfPointsForTriangleLumpedIntegration()) {
        const auto scheme = MakeLumpedIntegrationScheme(number);

        ExpectLocalCoordinatesAreInRange(scheme->GetIntegrationPoints(), 0.0, 1.0, std::size_t{0});
        ExpectLocalCoordinatesAreInRange(scheme->GetIntegrationPoints(), 0.0, 1.0, std::size_t{1});
        ExpectLocalCoordinatesAreZero(scheme->GetIntegrationPoints(), std::size_t{2});
        ExpectLocalCoordinatesIncludeRangeBounds(scheme->GetIntegrationPoints(), 0.0, 1.0, std::size_t{0});
        ExpectLocalCoordinatesIncludeRangeBounds(scheme->GetIntegrationPoints(), 0.0, 1.0, std::size_t{1});
    }
}

KRATOS_TEST_CASE_IN_SUITE(CorrectWeightsFromTriangle6LumpedSchemes, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const std::vector<double> expected_weights{0.0328638, 0.0328638, 0.0328638, 0.300469,
                                               0.300469,  0.300469}; // regression values see core test test_triangle_2d_6
    auto                      actual_weights = ComputeWeightsForLumpedIntegration(6);
    KRATOS_EXPECT_VECTOR_NEAR(expected_weights, actual_weights, 1.0E-6)
}

KRATOS_TEST_CASE_IN_SUITE(CorrectWeightsFromTriangle3LumpedSchemes, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const std::vector<double> expected_weights{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
    auto                      actual_weights = ComputeWeightsForLumpedIntegration(3);
    KRATOS_EXPECT_VECTOR_NEAR(expected_weights, actual_weights, 1.0E-6)
}

KRATOS_TEST_CASE_IN_SUITE(CorrectWeightsFromQuadrilateral4LumpedSchemes, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const std::vector<double> expected_weights{1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0};
    auto                      actual_weights = ComputeWeightsForLumpedIntegration(4);
    KRATOS_EXPECT_VECTOR_NEAR(expected_weights, actual_weights, 1.0E-6)
}

KRATOS_TEST_CASE_IN_SUITE(CorrectWeightsFromQuadrilateral8LumpedSchemes, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const std::vector<double> expected_weights{
        0.0394737, 0.0394737, 0.0394737, 0.0394737, 0.210526,
        0.210526,  0.210526,  0.210526}; // regression values see core test test_quadrilateral_2d_8
    auto actual_weights = ComputeWeightsForLumpedIntegration(8);
    KRATOS_EXPECT_VECTOR_NEAR(expected_weights, actual_weights, 1.0E-6)
}

} // namespace Kratos::Testing