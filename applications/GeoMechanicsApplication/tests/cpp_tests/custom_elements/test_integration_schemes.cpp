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
#include "tests/cpp_tests/test_utilities.h"

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

double SumOfWeights(const Geo::IntegrationPointVectorType& rIntegrationPoints)
{
    auto weights = std::vector<double>{};
    weights.reserve(rIntegrationPoints.size());
    std::ranges::transform(rIntegrationPoints, std::back_inserter(weights),
                           [](const auto& rIntegrationPoint) { return rIntegrationPoint.Weight(); });
    return std::accumulate(weights.cbegin(), weights.cend(), 0.0);
}

void ExpectLocalCoordinatesIncludeRangeBound(const Geo::IntegrationPointVectorType& rIntegrationPoints,
                                             double      Bound,
                                             std::size_t DirectionIndex)
{
    auto is_at_bound_of_iso_coord = [Bound, DirectionIndex](const auto& rPoint) {
        return std::abs(rPoint[DirectionIndex] - Bound) <= Testing::Defaults::absolute_tolerance;
    };
    KRATOS_EXPECT_TRUE(std::any_of(rIntegrationPoints.begin(), rIntegrationPoints.end(), is_at_bound_of_iso_coord))
}

void ExpectLocalCoordinatesAreInRange(const Geo::IntegrationPointVectorType& rIntegrationPoints,
                                      double                                 IsoLowerBound,
                                      double                                 IsoUpperBound,
                                      std::size_t                            DirectionIndex)
{
    auto iso_coord_is_in_range = [IsoLowerBound, IsoUpperBound, DirectionIndex](const auto& rPoint) {
        return rPoint[DirectionIndex] - IsoUpperBound <= Testing::Defaults::absolute_tolerance &&
               IsoLowerBound - rPoint[DirectionIndex] <= Testing::Defaults::absolute_tolerance;
    };
    KRATOS_EXPECT_TRUE(std::all_of(rIntegrationPoints.begin(), rIntegrationPoints.end(), iso_coord_is_in_range))
}

void ExpectLocalCoordinatesAreZero(const Geo::IntegrationPointVectorType& rIntegrationPoints, const size_t DirectionIndex)
{
    auto iso_coordinate_must_be_near_zero = [DirectionIndex](const auto& rPoint) {
        return (std::abs(rPoint[DirectionIndex]) <= Testing::Defaults::absolute_tolerance);
    };
    KRATOS_EXPECT_TRUE(std::all_of(rIntegrationPoints.begin(), rIntegrationPoints.end(), iso_coordinate_must_be_near_zero))
}

std::vector<double> ComputeWeightsForLumpedIntegration(const size_t NumberOfPoints)
{
    const auto          p_scheme           = MakeLumpedIntegrationScheme(NumberOfPoints);
    auto                integration_points = p_scheme->GetIntegrationPoints();
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
    const auto     p_scheme         = LobattoIntegrationScheme{number_of_points};

    KRATOS_EXPECT_NE(dynamic_cast<const IntegrationScheme*>(&p_scheme), nullptr);
}

class ParametrizedInterfaceIntegrationMethodSuite
    : public ::testing::TestWithParam<std::tuple<std::size_t, std::function<std::unique_ptr<IntegrationScheme>(std::size_t)>>>
{
};

TEST_P(ParametrizedInterfaceIntegrationMethodSuite, NumberOfIntegrationPointsMatchesTheNumberOfPointsGivenAtConstructionTime)
{
    const auto& [number, scheme_creator] = GetParam();
    const auto p_scheme                  = scheme_creator(number);

    KRATOS_EXPECT_EQ(p_scheme->GetNumberOfIntegrationPoints(), number);
    KRATOS_EXPECT_EQ(p_scheme->GetIntegrationPoints().size(), number);
}

INSTANTIATE_TEST_SUITE_P(KratosGeoMechanicsFastSuiteWithoutKernel,
                         ParametrizedInterfaceIntegrationMethodSuite,
                         ::testing::Values(std::make_tuple(std::size_t{2}, MakeLobattoIntegrationScheme),
                                           std::make_tuple(std::size_t{3}, MakeLobattoIntegrationScheme),
                                           std::make_tuple(std::size_t{3}, MakeLumpedIntegrationScheme),
                                           std::make_tuple(std::size_t{6}, MakeLumpedIntegrationScheme),
                                           std::make_tuple(std::size_t{4}, MakeLumpedIntegrationScheme),
                                           std::make_tuple(std::size_t{8}, MakeLumpedIntegrationScheme)));

KRATOS_TEST_CASE_IN_SUITE(SumOfIntegrationPointWeightsOfAllSupportedLobattoSchemesEqualsTwo,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    for (auto number : SupportedNumbersOfPointsForLobattoIntegration()) {
        const auto p_scheme = MakeLobattoIntegrationScheme(number);
        KRATOS_EXPECT_RELATIVE_NEAR(SumOfWeights(p_scheme->GetIntegrationPoints()), 2.0, Defaults::relative_tolerance)
    }
}

KRATOS_TEST_CASE_IN_SUITE(PointsOfAllSupportedLobattoSchemesMustBeInRangeAndIncludeBounds,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    for (auto number : SupportedNumbersOfPointsForLobattoIntegration()) {
        const auto scheme = MakeLobattoIntegrationScheme(number);

        ExpectLocalCoordinatesAreInRange(scheme->GetIntegrationPoints(), -1.0, 1.0, std::size_t{0});
        ExpectLocalCoordinatesIncludeRangeBound(scheme->GetIntegrationPoints(), -1.0, std::size_t{0});
        ExpectLocalCoordinatesIncludeRangeBound(scheme->GetIntegrationPoints(), 1.0, std::size_t{0});
        ExpectLocalCoordinatesAreZero(scheme->GetIntegrationPoints(), std::size_t{1});
        ExpectLocalCoordinatesAreZero(scheme->GetIntegrationPoints(), std::size_t{2});
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
    const auto     p_scheme         = LumpedIntegrationScheme{number_of_points};

    KRATOS_EXPECT_NE(dynamic_cast<const IntegrationScheme*>(&p_scheme), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(SumOfIntegrationPointWeightsOfAllSupportedLumpedSchemesEqualsArea,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    constexpr auto area_quadrilateral = 4.0;
    for (auto number : SupportedNumbersOfPointsForQuadrilateralLumpedIntegration()) {
        const auto p_scheme = MakeLumpedIntegrationScheme(number);
        KRATOS_EXPECT_RELATIVE_NEAR(SumOfWeights(p_scheme->GetIntegrationPoints()),
                                    area_quadrilateral, Defaults::relative_tolerance)
    }

    constexpr auto area_triangle = 0.5;
    for (auto number : SupportedNumbersOfPointsForTriangleLumpedIntegration()) {
        const auto p_scheme = MakeLumpedIntegrationScheme(number);
        KRATOS_EXPECT_RELATIVE_NEAR(SumOfWeights(p_scheme->GetIntegrationPoints()), area_triangle,
                                    Defaults::relative_tolerance)
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
        const auto p_scheme = MakeLumpedIntegrationScheme(number);

        ExpectLocalCoordinatesAreInRange(p_scheme->GetIntegrationPoints(), -1.0, 1.0, std::size_t{0});
        ExpectLocalCoordinatesAreInRange(p_scheme->GetIntegrationPoints(), -1.0, 1.0, std::size_t{1});
        ExpectLocalCoordinatesIncludeRangeBound(p_scheme->GetIntegrationPoints(), -1.0, std::size_t{0});
        ExpectLocalCoordinatesIncludeRangeBound(p_scheme->GetIntegrationPoints(), 1.0, std::size_t{0});
        ExpectLocalCoordinatesIncludeRangeBound(p_scheme->GetIntegrationPoints(), -1.0, std::size_t{1});
        ExpectLocalCoordinatesIncludeRangeBound(p_scheme->GetIntegrationPoints(), 1.0, std::size_t{1});
        ExpectLocalCoordinatesAreZero(p_scheme->GetIntegrationPoints(), std::size_t{2});
    }
}

KRATOS_TEST_CASE_IN_SUITE(PointsOfAllSupportedTriangleLumpedSchemesMustBeInRangeAndIncludeBounds,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    for (auto number : SupportedNumbersOfPointsForTriangleLumpedIntegration()) {
        const auto p_scheme = MakeLumpedIntegrationScheme(number);

        ExpectLocalCoordinatesAreInRange(p_scheme->GetIntegrationPoints(), 0.0, 1.0, std::size_t{0});
        ExpectLocalCoordinatesAreInRange(p_scheme->GetIntegrationPoints(), 0.0, 1.0, std::size_t{1});
        ExpectLocalCoordinatesIncludeRangeBound(p_scheme->GetIntegrationPoints(), 0.0, std::size_t{0});
        ExpectLocalCoordinatesIncludeRangeBound(p_scheme->GetIntegrationPoints(), 1.0, std::size_t{0});
        ExpectLocalCoordinatesIncludeRangeBound(p_scheme->GetIntegrationPoints(), 0.0, std::size_t{1});
        ExpectLocalCoordinatesIncludeRangeBound(p_scheme->GetIntegrationPoints(), 1.0, std::size_t{1});
        ExpectLocalCoordinatesAreZero(p_scheme->GetIntegrationPoints(), std::size_t{2});
    }
}

KRATOS_TEST_CASE_IN_SUITE(CorrectWeightsFromTriangle6LumpedSchemes, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    std::vector<double> expected_weights{0.0328638, 0.0328638, 0.0328638, 0.300469,
                                         0.300469,  0.300469}; // regression values see core test test_triangle_2d_6
    std::transform(expected_weights.begin(), expected_weights.end(), expected_weights.begin(),
                   [](const auto& rWeight) {
        return rWeight * 0.5; // Adjust weights for the area of the triangle
    });
    const auto     actual_weights = ComputeWeightsForLumpedIntegration(6);
    constexpr auto tolerance      = 1e-5;
    KRATOS_EXPECT_VECTOR_NEAR(expected_weights, actual_weights, tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(CorrectWeightsFromTriangle3LumpedSchemes, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    std::vector<double> expected_weights{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
    std::transform(expected_weights.begin(), expected_weights.end(), expected_weights.begin(),
                   [](const auto& rWeight) {
        return rWeight * 0.5; // Adjust weights for the area of the triangle
    });
    const auto     actual_weights = ComputeWeightsForLumpedIntegration(3);
    constexpr auto tolerance      = 1e-5;
    KRATOS_EXPECT_VECTOR_NEAR(expected_weights, actual_weights, tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(CorrectWeightsFromQuadrilateral4LumpedSchemes, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    std::vector<double> expected_weights{1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0};
    std::transform(expected_weights.begin(), expected_weights.end(), expected_weights.begin(),
                   [](const auto& rWeight) {
        return rWeight * 4.0; // Adjust weights for the area of the quadrilateral
    });
    const auto     actual_weights = ComputeWeightsForLumpedIntegration(4);
    constexpr auto tolerance      = 1e-5;
    KRATOS_EXPECT_VECTOR_NEAR(expected_weights, actual_weights, tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(CorrectWeightsFromQuadrilateral8LumpedSchemes, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    std::vector<double> expected_weights{
        0.0394737, 0.0394737, 0.0394737, 0.0394737, 0.210526,
        0.210526,  0.210526,  0.210526}; // regression lumping values see core test test_quadrilateral_2d_8
    std::transform(expected_weights.begin(), expected_weights.end(), expected_weights.begin(),
                   [](const auto& rWeight) {
        return rWeight * 4.0; // Adjust weights for the area of the quadrilateral
    });
    const auto     actual_weights = ComputeWeightsForLumpedIntegration(8);
    constexpr auto tolerance      = 1e-5;
    KRATOS_EXPECT_VECTOR_NEAR(expected_weights, actual_weights, tolerance)
}

} // namespace Kratos::Testing