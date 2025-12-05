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
#include "test_setup_utilities/element_setup_utilities.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

#include <boost/numeric/ublas/assignment.hpp>
#include <numbers>

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(NodalExtrapolator_GivesCorrectExtrapolationMatrix_For2D2NLine,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto nodes = ElementSetupUtilities::GenerateNodes({{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}});
    auto p_element = ElementSetupUtilities::Create2D2NElement(nodes, std::make_shared<Properties>());

    const LinearNodalExtrapolator nodal_extrapolator;
    const auto extrapolation_matrix = nodal_extrapolator.CalculateElementExtrapolationMatrix(*p_element);

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
    auto nodes = ElementSetupUtilities::GenerateNodes({{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.5, 0.0, 0.0}});
    auto p_element = ElementSetupUtilities::Create2D3NLineElement(nodes, std::make_shared<Properties>());

    const LinearNodalExtrapolator nodal_extrapolator;
    const auto extrapolation_matrix = nodal_extrapolator.CalculateElementExtrapolationMatrix(*p_element);

    // clang-format off
    const auto expected_extrapolation_matrix = UblasUtilities::CreateMatrix(
        {{0.5 * (1 + std::numbers::sqrt3), 0.5 * (1 - std::numbers::sqrt3)},
         {0.5 * (1 - std::numbers::sqrt3), 0.5 * (1 + std::numbers::sqrt3)},
         {0.5, 0.5}});
    // clang-format on

    KRATOS_EXPECT_MATRIX_NEAR(extrapolation_matrix, expected_extrapolation_matrix, Defaults::absolute_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(NodalExtrapolator_GivesCorrectExtrapolationMatrix_For2Plus2LineInterface,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto nodes = ElementSetupUtilities::GenerateNodes(
        {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 0.0, 0.5}, {1.0, 0.0, 0.5}});
    auto p_element =
        ElementSetupUtilities::Create2D4NInterfaceElement(nodes, std::make_shared<Properties>());

    const LinearNodalExtrapolator nodal_extrapolator;
    auto extrapolation_matrix = nodal_extrapolator.CalculateElementExtrapolationMatrix(*p_element);

    // clang-format off
    const auto expected_extrapolation_matrix = UblasUtilities::CreateMatrix({{ 1.0, 0.0},
                                                                             {0.0, 1.0},
                                                                             {1.0, 0.0},
                                                                             {0.0, 1.0}});
    // clang-format on

    KRATOS_EXPECT_MATRIX_NEAR(extrapolation_matrix, expected_extrapolation_matrix, Defaults::absolute_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(NodalExtrapolator_GivesCorrectExtrapolationMatrix_For2D3NTriangle,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto nodes = ElementSetupUtilities::GenerateNodes({{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}});
    auto p_element = ElementSetupUtilities::Create2D3NElement(nodes, std::make_shared<Properties>());

    const LinearNodalExtrapolator nodal_extrapolator;
    auto extrapolation_matrix = nodal_extrapolator.CalculateElementExtrapolationMatrix(*p_element);

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
    auto nodes = ElementSetupUtilities::GenerateNodes(
        {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.5, 0.0, 0.0}, {0.5, 0.5, 0.0}, {0.0, 0.5, 0.0}});
    auto p_element = ElementSetupUtilities::Create2D6NElement(nodes, std::make_shared<Properties>());

    const LinearNodalExtrapolator nodal_extrapolator;
    auto extrapolation_matrix = nodal_extrapolator.CalculateElementExtrapolationMatrix(*p_element);

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

KRATOS_TEST_CASE_IN_SUITE(NodalExtrapolator_GivesCorrectExtrapolationMatrix_For3Plus3TriangularInterface,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto nodes = ElementSetupUtilities::GenerateNodes(
        {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 0.1}, {1.0, 0.0, 0.1}, {0.0, 1.0, 0.1}});
    auto p_element =
        ElementSetupUtilities::Create3D6NInterfaceElement(nodes, std::make_shared<Properties>());

    const LinearNodalExtrapolator nodal_extrapolator;
    auto extrapolation_matrix = nodal_extrapolator.CalculateElementExtrapolationMatrix(*p_element);

    // clang-format off
    const auto expected_extrapolation_matrix = UblasUtilities::CreateMatrix({{1.0, 0.0, 0.0},
                                                                            {0.0, 1.0, 0.0},
                                                                            {0.0, 0.0, 1.0},
                                                                            {1.0, 0.0, 0.0},
                                                                            {0.0, 1.0, 0.0},
                                                                            {0.0, 0.0, 1.0}});
    // clang-format on

    KRATOS_EXPECT_MATRIX_NEAR(extrapolation_matrix, expected_extrapolation_matrix, Defaults::absolute_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(NodalExtrapolator_GivesCorrectExtrapolationMatrix_For4Plus4QuadrilateralInterface,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto nodes = ElementSetupUtilities::GenerateNodes({{0.0, 0.0, 0.0},
                                                       {1.0, 0.0, 0.0},
                                                       {1.0, 1.0, 0.0},
                                                       {0.0, 1.0, 0.0},
                                                       {0.0, 0.0, 0.1},
                                                       {1.0, 0.0, 0.1},
                                                       {1.0, 1.0, 0.1},
                                                       {0.0, 1.0, 0.1}});
    auto p_element =
        ElementSetupUtilities::Create3D8NInterfaceElement(nodes, std::make_shared<Properties>(), 1);

    const LinearNodalExtrapolator nodal_extrapolator;
    auto extrapolation_matrix = nodal_extrapolator.CalculateElementExtrapolationMatrix(*p_element);

    // clang-format off
    const auto expected_extrapolation_matrix = UblasUtilities::CreateMatrix({{1.0, 0.0, 0.0, 0.0},
                                                                            {0.0, 1.0, 0.0, 0.0},
                                                                            {0.0, 0.0, 1.0, 0.0},
                                                                            {0.0, 0.0, 0.0, 1.0},
                                                                            {1.0, 0.0, 0.0, 0.0},
                                                                            {0.0, 1.0, 0.0, 0.0},
                                                                            {0.0, 0.0, 1.0, 0.0},
                                                                            {0.0, 0.0, 0.0, 1.0}});
    // clang-format on

    KRATOS_EXPECT_MATRIX_NEAR(extrapolation_matrix, expected_extrapolation_matrix, Defaults::absolute_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(NodalExtrapolator_GivesCorrectExtrapolationMatrix_For2D4NQuadrilateral,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto nodes = ElementSetupUtilities::GenerateNodes(
        {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0}});
    auto p_element = ElementSetupUtilities::Create2D4NElement(nodes, std::make_shared<Properties>());

    const LinearNodalExtrapolator nodal_extrapolator;
    auto extrapolation_matrix = nodal_extrapolator.CalculateElementExtrapolationMatrix(*p_element);

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
    auto nodes = ElementSetupUtilities::GenerateNodes({{0.0, 0.0, 0.0},
                                                       {1.0, 0.0, 0.0},
                                                       {1.0, 1.0, 0.0},
                                                       {0.0, 1.0, 0.0},
                                                       {0.5, 0.0, 0.0},
                                                       {1.0, 0.5, 0.0},
                                                       {0.5, 1.0, 0.0},
                                                       {0.0, 0.5, 0.0}});
    auto p_element = ElementSetupUtilities::Create2D8NElement(nodes, std::make_shared<Properties>());

    const LinearNodalExtrapolator nodal_extrapolator;
    auto extrapolation_matrix = nodal_extrapolator.CalculateElementExtrapolationMatrix(*p_element);

    // clang-format off
    const auto expected_extrapolation_matrix = UblasUtilities::CreateMatrix({{ 1.866025, -0.5,       0.133974, -0.5     },
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
    auto nodes = ElementSetupUtilities::GenerateNodes(
        {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}});
    auto p_element = ElementSetupUtilities::Create3D4NElement(nodes, std::make_shared<Properties>());

    const LinearNodalExtrapolator nodal_extrapolator;
    auto extrapolation_matrix = nodal_extrapolator.CalculateElementExtrapolationMatrix(*p_element);

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
    auto nodes = ElementSetupUtilities::GenerateNodes({{0.0, 0.0, 0.0},
                                                       {1.0, 0.0, 0.0},
                                                       {0.0, 1.0, 0.0},
                                                       {0.0, 0.0, 1.0},
                                                       {0.5, 0.0, 0.0},
                                                       {0.5, 0.5, 0.0},
                                                       {0.0, 0.5, 0.0},
                                                       {0.0, 0.0, 0.5},
                                                       {0.5, 0.0, 0.5},
                                                       {0.0, 0.5, 0.5}});
    auto p_element = ElementSetupUtilities::Create3D10NElement(nodes, std::make_shared<Properties>());

    const LinearNodalExtrapolator nodal_extrapolator;
    auto extrapolation_matrix = nodal_extrapolator.CalculateElementExtrapolationMatrix(*p_element);

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
    auto nodes = ElementSetupUtilities::GenerateNodes({{0.0, 0.0, 0.0},
                                                       {1.0, 0.0, 0.0},
                                                       {1.0, 1.0, 0.0},
                                                       {0.0, 1.0, 0.0},
                                                       {0.0, 0.0, 1.0},
                                                       {1.0, 0.0, 1.0},
                                                       {1.0, 1.0, 1.0},
                                                       {0.0, 1.0, 1.0}});
    auto p_element = ElementSetupUtilities::Create3D8NElement(nodes, std::make_shared<Properties>());

    const LinearNodalExtrapolator nodal_extrapolator;
    auto extrapolation_matrix = nodal_extrapolator.CalculateElementExtrapolationMatrix(*p_element);

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
    auto nodes = ElementSetupUtilities::GenerateNodes(
        {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0},
         {1.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0}, {0.5, 0.0, 0.0}, {1.0, 0.5, 0.0},
         {0.5, 1.0, 0.0}, {0.0, 0.5, 0.0}, {0.0, 0.0, 0.5}, {1.0, 0.0, 0.5}, {1.0, 1.0, 0.5},
         {0.0, 1.0, 0.5}, {0.5, 0.0, 1.0}, {1.0, 0.5, 1.0}, {0.5, 1.0, 1.0}, {0.0, 0.5, 1.0}});
    auto p_element = ElementSetupUtilities::Create3D20NElement(nodes, std::make_shared<Properties>());

    const LinearNodalExtrapolator nodal_extrapolator;
    auto extrapolation_matrix = nodal_extrapolator.CalculateElementExtrapolationMatrix(*p_element);

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
