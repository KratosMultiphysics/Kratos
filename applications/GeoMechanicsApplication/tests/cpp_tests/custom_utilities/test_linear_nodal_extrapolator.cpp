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
//                   Aron Noordam
//

#include "custom_utilities/linear_nodal_extrapolator.h"
#include "custom_utilities/ublas_utilities.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_2d_3.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"
#include "geometries/tetrahedra_3d_10.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"
#include "test_setup_utilities/element_setup_utilities.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

#include <boost/numeric/ublas/assignment.hpp>
#include <numbers>

namespace Kratos::Testing
{

using namespace Kratos;

KRATOS_TEST_CASE_IN_SUITE(NodalExtrapolator_GivesCorrectExtrapolationMatrix_For2D2NLine,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto nodes    = ElementSetupUtilities::GenerateNodes({{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}});
    const auto geometry = Line2D2<Node>{nodes};

    const auto     nodal_extrapolator   = LinearNodalExtrapolator{};
    constexpr auto integration_method   = GeometryData::IntegrationMethod::GI_GAUSS_2;
    const auto     extrapolation_matrix = nodal_extrapolator.CalculateElementExtrapolationMatrix(
        geometry, geometry.IntegrationPoints(integration_method));

    // clang-format off
    const auto expected_extrapolation_matrix = UblasUtilities::CreateMatrix(
        {{0.5 * (1 + std::numbers::sqrt3), 0.5 * (1 - std::numbers::sqrt3)},
         {0.5 * (1 - std::numbers::sqrt3), 0.5 * (1 + std::numbers::sqrt3)}});
    // clang-format on

    KRATOS_EXPECT_MATRIX_NEAR(extrapolation_matrix, expected_extrapolation_matrix, Defaults::absolute_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(NodalExtrapolator_GivesCorrectExtrapolationMatrix_For2D3NLine,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto nodes =
        ElementSetupUtilities::GenerateNodes({{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.5, 0.0, 0.0}});
    const auto geometry = Line2D3<Node>{nodes};

    const auto     nodal_extrapolator   = LinearNodalExtrapolator{};
    constexpr auto integration_method   = GeometryData::IntegrationMethod::GI_GAUSS_2;
    const auto     extrapolation_matrix = nodal_extrapolator.CalculateElementExtrapolationMatrix(
        geometry, geometry.IntegrationPoints(integration_method));

    // clang-format off
    const auto expected_extrapolation_matrix = UblasUtilities::CreateMatrix(
        {{0.5 * (1 + std::numbers::sqrt3), 0.5 * (1 - std::numbers::sqrt3)},
         {0.5 * (1 - std::numbers::sqrt3), 0.5 * (1 + std::numbers::sqrt3)},
         {0.5, 0.5}});
    // clang-format on

    KRATOS_EXPECT_MATRIX_NEAR(extrapolation_matrix, expected_extrapolation_matrix, Defaults::absolute_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(NodalExtrapolator_GivesCorrectExtrapolationMatrix_For2D3NTriangle,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto nodes =
        ElementSetupUtilities::GenerateNodes({{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}});
    const auto geometry = Triangle2D3<Node>{nodes};

    const auto     nodal_extrapolator   = LinearNodalExtrapolator{};
    constexpr auto integration_method   = GeometryData::IntegrationMethod::GI_GAUSS_2;
    const auto     extrapolation_matrix = nodal_extrapolator.CalculateElementExtrapolationMatrix(
        geometry, geometry.IntegrationPoints(integration_method));

    // clang-format off
    const auto expected_extrapolation_matrix = UblasUtilities::CreateMatrix({{ 5.0/3.0, -1.0/3.0, -1.0/3.0},
                                                                             {-1.0/3.0,  5.0/3.0, -1.0/3.0},
                                                                             {-1.0/3.0, -1.0/3.0,  5.0/3.0}});
    // clang-format on

    KRATOS_EXPECT_MATRIX_NEAR(extrapolation_matrix, expected_extrapolation_matrix, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(NodalExtrapolator_GivesCorrectExtrapolationMatrix_For2D6NTriangle,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto nodes = ElementSetupUtilities::GenerateNodes(
        {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.5, 0.0, 0.0}, {0.5, 0.5, 0.0}, {0.0, 0.5, 0.0}});
    const auto geometry = Triangle2D6<Node>{nodes};

    const auto     nodal_extrapolator   = LinearNodalExtrapolator{};
    constexpr auto integration_method   = GeometryData::IntegrationMethod::GI_GAUSS_2;
    const auto     extrapolation_matrix = nodal_extrapolator.CalculateElementExtrapolationMatrix(
        geometry, geometry.IntegrationPoints(integration_method));

    // clang-format off
    const auto expected_extrapolation_matrix = UblasUtilities::CreateMatrix({{ 5.0/3.0, -1.0/3.0, -1.0/3.0},
                                                                             {-1.0/3.0,  5.0/3.0, -1.0/3.0},
                                                                             {-1.0/3.0, -1.0/3.0,  5.0/3.0},
                                                                             { 2.0/3.0,  2.0/3.0, -1.0/3.0},
                                                                             {-1.0/3.0,  2.0/3.0,  2.0/3.0},
                                                                             { 2.0/3.0, -1.0/3.0,  2.0/3.0}});
    // clang-format on

    KRATOS_EXPECT_MATRIX_NEAR(extrapolation_matrix, expected_extrapolation_matrix, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(NodalExtrapolator_GivesCorrectExtrapolationMatrix_For2D4NQuadrilateral,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto nodes = ElementSetupUtilities::GenerateNodes(
        {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0}});
    const auto geometry = Quadrilateral2D4<Node>{nodes};

    const auto     nodal_extrapolator   = LinearNodalExtrapolator{};
    constexpr auto integration_method   = GeometryData::IntegrationMethod::GI_GAUSS_2;
    const auto     extrapolation_matrix = nodal_extrapolator.CalculateElementExtrapolationMatrix(
        geometry, geometry.IntegrationPoints(integration_method));

    // clang-format off
    const auto expected_extrapolation_matrix = UblasUtilities::CreateMatrix({{ 1.866025, -0.5,       0.133974, -0.5     },
                                                                             {-0.5,       1.866025, -0.5,       0.133974},
                                                                             { 0.133974, -0.5,       1.866025, -0.5     },
                                                                             {-0.5,       0.133974, -0.5,       1.866025}});
    // clang-format on

    KRATOS_EXPECT_MATRIX_NEAR(extrapolation_matrix, expected_extrapolation_matrix, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(NodalExtrapolator_GivesCorrectExtrapolationMatrix_For2D8NQuadrilateral,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto nodes    = ElementSetupUtilities::GenerateNodes({{0.0, 0.0, 0.0},
                                                                {1.0, 0.0, 0.0},
                                                                {1.0, 1.0, 0.0},
                                                                {0.0, 1.0, 0.0},
                                                                {0.5, 0.0, 0.0},
                                                                {1.0, 0.5, 0.0},
                                                                {0.5, 1.0, 0.0},
                                                                {0.0, 0.5, 0.0}});
    const auto geometry = Quadrilateral2D8<Node>{nodes};

    const auto     nodal_extrapolator   = LinearNodalExtrapolator{};
    constexpr auto integration_method   = GeometryData::IntegrationMethod::GI_GAUSS_2;
    const auto     extrapolation_matrix = nodal_extrapolator.CalculateElementExtrapolationMatrix(
        geometry, geometry.IntegrationPoints(integration_method));

    // clang-format off
    const auto expected_extrapolation_matrix = UblasUtilities::CreateMatrix(
                                    {{ 1.866025, -0.5,      0.133974, -0.5     },
                                     {-0.5,       1.866025, -0.5,       0.133974},
                                     { 0.133974, -0.5,       1.866025, -0.5     },
                                     {-0.5,       0.133974, -0.5,       1.866025},
                                     { 0.683013,  0.683013, -0.183013, -0.183013},
                                     {-0.183013,  0.683013,  0.683013, -0.183013},
                                     {-0.183013, -0.183013,  0.683013,  0.683013},
                                     { 0.683013, -0.183013, -0.183013,  0.683013}});
    // clang-format on

    KRATOS_EXPECT_MATRIX_NEAR(extrapolation_matrix, expected_extrapolation_matrix, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(NodalExtrapolator_GivesCorrectExtrapolationMatrix_For3D4NKratos_Tetrahedra,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto nodes = ElementSetupUtilities::GenerateNodes(
        {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}});
    const auto geometry = Tetrahedra3D4<Node>{nodes};

    const auto     nodal_extrapolator   = LinearNodalExtrapolator{};
    constexpr auto integration_method   = GeometryData::IntegrationMethod::GI_GAUSS_2;
    const auto     extrapolation_matrix = nodal_extrapolator.CalculateElementExtrapolationMatrix(
        geometry, geometry.IntegrationPoints(integration_method));

    // clang-format off
    const auto expected_extrapolation_matrix = UblasUtilities::CreateMatrix({{-0.309017, -0.309017, -0.309017,  1.927051},
                                                                             { 1.927051, -0.309017, -0.309017, -0.309017},
                                                                             {-0.309017,  1.927051, -0.309017, -0.309017},
                                                                             {-0.309017, -0.309017,  1.927051, -0.309017}});
    // clang-format on

    KRATOS_EXPECT_MATRIX_NEAR(extrapolation_matrix, expected_extrapolation_matrix, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(NodalExtrapolator_GivesCorrectExtrapolationMatrix_For3D10NKratos_Tetrahedra,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto nodes    = ElementSetupUtilities::GenerateNodes({{0.0, 0.0, 0.0},
                                                                {1.0, 0.0, 0.0},
                                                                {0.0, 1.0, 0.0},
                                                                {0.0, 0.0, 1.0},
                                                                {0.5, 0.0, 0.0},
                                                                {0.5, 0.5, 0.0},
                                                                {0.0, 0.5, 0.0},
                                                                {0.0, 0.0, 0.5},
                                                                {0.5, 0.0, 0.5},
                                                                {0.0, 0.5, 0.5}});
    const auto geometry = Tetrahedra3D10<Node>{nodes};

    const auto     nodal_extrapolator   = LinearNodalExtrapolator{};
    constexpr auto integration_method   = GeometryData::IntegrationMethod::GI_GAUSS_2;
    const auto     extrapolation_matrix = nodal_extrapolator.CalculateElementExtrapolationMatrix(
        geometry, geometry.IntegrationPoints(integration_method));

    // clang-format off
    const auto expected_extrapolation_matrix = UblasUtilities::CreateMatrix({{-0.309017, -0.309017, -0.309017,  1.927051},
                                                                             { 1.927051, -0.309017, -0.309017, -0.309017},
                                                                             {-0.309017,  1.927051, -0.309017, -0.309017},
                                                                             {-0.309017, -0.309017,  1.927051, -0.309017},

                                                                             { 0.809017, -0.309017, -0.309017,  0.809017},
                                                                             { 0.809017,  0.809017, -0.309017, -0.309017},
                                                                             {-0.309017,  0.809017, -0.309017,  0.809017},
                                                                             {-0.309017, -0.309017,  0.809017,  0.809017},
                                                                             { 0.809017, -0.309017,  0.809017, -0.309017},
                                                                             {-0.309017,  0.809017,  0.809017, -0.309017}});
    // clang-format on

    KRATOS_EXPECT_MATRIX_NEAR(extrapolation_matrix, expected_extrapolation_matrix, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(NodalExtrapolator_GivesCorrectExtrapolationMatrix_For3D8NKratos_Hexahedra,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto nodes    = ElementSetupUtilities::GenerateNodes({{0.0, 0.0, 0.0},
                                                                {1.0, 0.0, 0.0},
                                                                {1.0, 1.0, 0.0},
                                                                {0.0, 1.0, 0.0},
                                                                {0.0, 0.0, 1.0},
                                                                {1.0, 0.0, 1.0},
                                                                {1.0, 1.0, 1.0},
                                                                {0.0, 1.0, 1.0}});
    const auto geometry = Hexahedra3D8<Node>{nodes};

    const auto     nodal_extrapolator   = LinearNodalExtrapolator{};
    constexpr auto integration_method   = GeometryData::IntegrationMethod::GI_GAUSS_2;
    const auto     extrapolation_matrix = nodal_extrapolator.CalculateElementExtrapolationMatrix(
        geometry, geometry.IntegrationPoints(integration_method));

    // clang-format off
    const auto expected_extrapolation_matrix = UblasUtilities::CreateMatrix({{ 2.549038, -0.683013,  0.183013, -0.683013, -0.683013,  0.183013, -0.049038,  0.183013},
                                                                             {-0.683013,  2.549038, -0.683013,  0.183013,  0.183013, -0.683013,  0.183013, -0.049038},
                                                                             { 0.183013, -0.683013,  2.549038, -0.683013, -0.049038,  0.183013, -0.683013,  0.183013},
                                                                             {-0.683013,  0.183013, -0.683013,  2.549038,  0.183013, -0.049038,  0.183013, -0.683013},
                                                                             {-0.683013,  0.183013, -0.049038,  0.183013,  2.549038, -0.683013,  0.183013, -0.683013},
                                                                             { 0.183013, -0.683013,  0.183013, -0.049038, -0.683013,  2.549038, -0.683013,  0.183013},
                                                                             {-0.049038,  0.183013, -0.683013,  0.183013,  0.183013, -0.683013,  2.549038, -0.683013},
                                                                             { 0.183013, -0.049038,  0.183013, -0.683013, -0.683013,  0.183013, -0.683013,  2.549038}});
    // clang-format on

    KRATOS_EXPECT_MATRIX_NEAR(extrapolation_matrix, expected_extrapolation_matrix, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(NodalExtrapolator_GivesCorrectExtrapolationMatrix_For3D20NKratos_Hexahedra,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto nodes = ElementSetupUtilities::GenerateNodes(
        {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0},
         {1.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0}, {0.5, 0.0, 0.0}, {1.0, 0.5, 0.0},
         {0.5, 1.0, 0.0}, {0.0, 0.5, 0.0}, {0.0, 0.0, 0.5}, {1.0, 0.0, 0.5}, {1.0, 1.0, 0.5},
         {0.0, 1.0, 0.5}, {0.5, 0.0, 1.0}, {1.0, 0.5, 1.0}, {0.5, 1.0, 1.0}, {0.0, 0.5, 1.0}});
    const auto geometry = Hexahedra3D20<Node>{nodes};

    const auto     nodal_extrapolator   = LinearNodalExtrapolator{};
    constexpr auto integration_method   = GeometryData::IntegrationMethod::GI_GAUSS_2;
    const auto     extrapolation_matrix = nodal_extrapolator.CalculateElementExtrapolationMatrix(
        geometry, geometry.IntegrationPoints(integration_method));

    // clang-format off
    const auto expected_extrapolation_matrix = UblasUtilities::CreateMatrix({{ 2.549038,  -0.683013,   0.183013,  -0.683013,  -0.683013,   0.183013,  -0.049038,   0.183013},
                                                                             {-0.683013,   2.549038,  -0.683013,   0.183013,   0.183013,  -0.683013,   0.183013,  -0.049038},
                                                                             { 0.183013,  -0.683013,   2.549038,  -0.683013,  -0.049038,   0.183013,  -0.683013,   0.183013},
                                                                             {-0.683013,   0.183013,  -0.683013,   2.549038,   0.183013,  -0.049038,   0.183013,  -0.683013},
                                                                             {-0.683013,   0.183013,  -0.049038,   0.183013,   2.549038,  -0.683013,   0.183013,  -0.683013},
                                                                             { 0.183013,  -0.683013,   0.183013,  -0.049038,  -0.683013,   2.549038,  -0.683013,   0.183013},
                                                                             {-0.049038,   0.183013,  -0.683013,   0.183013,   0.183013,  -0.683013,   2.549038,  -0.683013},
                                                                             { 0.183013,  -0.049038,   0.183013,  -0.683013,  -0.683013,   0.183013,  -0.683013,   2.549038},

                                                                             { 0.9330125,  0.9330125, -0.25,      -0.25,      -0.25,      -0.25,       0.0669875,  0.0669875},
                                                                             {-0.25,       0.9330125,  0.9330125, -0.25,       0.0669875, -0.25,      -0.25,       0.0669875},
                                                                             {-0.25,      -0.25,       0.9330125,  0.9330125,  0.0669875,  0.0669875, -0.25,      -0.25     },
                                                                             { 0.9330125, -0.25,      -0.25,       0.9330125, -0.25,       0.0669875,  0.0669875, -0.25     },
                                                                             { 0.9330125, -0.25,       0.0669875, -0.25,       0.9330125, -0.25,       0.0669875, -0.25     },
                                                                             {-0.25,       0.9330125, -0.25,       0.0669875, -0.25,       0.9330125, -0.25,       0.0669875},
                                                                             { 0.0669875, -0.25,       0.9330125, -0.25,       0.0669875, -0.25,       0.9330125, -0.25     },
                                                                             {-0.25,       0.0669875, -0.25,       0.9330125, -0.25,       0.0669875, -0.25,       0.9330125},
                                                                             {-0.25,      -0.25,       0.0669875,  0.0669875,  0.9330125,  0.9330125, -0.25,      -0.25     },
                                                                             { 0.0669875, -0.25,      -0.25,       0.0669875, -0.25,       0.9330125,  0.9330125, -0.25     },
                                                                             { 0.0669875,  0.0669875, -0.25,      -0.25,      -0.25,      -0.25,       0.9330125,  0.9330125},
                                                                             {-0.25,       0.0669875,  0.0669875, -0.25,       0.9330125, -0.25,      -0.25,       0.9330125}});

    // clang-format on

    KRATOS_EXPECT_MATRIX_NEAR(extrapolation_matrix, expected_extrapolation_matrix, 1e-6)
}

} // namespace Kratos::Testing
