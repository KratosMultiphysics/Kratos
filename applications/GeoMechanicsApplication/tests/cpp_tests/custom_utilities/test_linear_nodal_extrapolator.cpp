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
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include <boost/numeric/ublas/assignment.hpp>

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(NodalExtrapolator_GivesCorrectExtrapolationMatrix_For2D2NLine,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Line2D2<Node> geometry(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0),
                           Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0));

    const LinearNodalExtrapolator nodal_extrapolator;

    // Used the integration method found in the 2d3 Timoshenko Beam element
    constexpr auto integration_method = GeometryData::IntegrationMethod::GI_GAUSS_3;
    auto           extrapolation_matrix =
        nodal_extrapolator.CalculateElementExtrapolationMatrix(geometry, integration_method);

    // clang-format off
    Matrix expected_extrapolation_matrix = ZeroMatrix(2, 3);
    expected_extrapolation_matrix <<= UblasUtilities::CreateMatrix({{ 0.923275, 0.444444,-0.367719},
                                                                    {-0.367719, 0.444444, 0.923275}});
    // clang-format on

    KRATOS_EXPECT_MATRIX_NEAR(extrapolation_matrix, expected_extrapolation_matrix, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(NodalExtrapolator_GivesCorrectExtrapolationMatrix_For2D3NLine,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Line2D3<Node> geometry(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0),
                           Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0),
                           Kratos::make_intrusive<Node>(3, 0.5, 0.0, 0.0));

    const LinearNodalExtrapolator nodal_extrapolator;

    // Used the integration method found in the 2d2 Timoshenko Beam element
    constexpr auto integration_method = GeometryData::IntegrationMethod::GI_GAUSS_5;
    auto           extrapolation_matrix =
        nodal_extrapolator.CalculateElementExtrapolationMatrix(geometry, integration_method);

    // clang-format off
    Matrix expected_extrapolation_matrix = ZeroMatrix(3, 5);
    expected_extrapolation_matrix <<= UblasUtilities::CreateMatrix(
        {{0.873934, 0.548036,-0.426666,-0.225144, 0.229839},
         {0.229839,-0.225144,-0.426666, 0.548036, 0.873934},
         {-0.098248,0.278248, 0.640000, 0.278248,-0.098248}});
    // clang-format on

    KRATOS_EXPECT_MATRIX_NEAR(extrapolation_matrix, expected_extrapolation_matrix, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(NodalExtrapolator_GivesCorrectExtrapolationMatrix_For2D3NTriangle,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Kratos::Triangle2D3<Node> geometry(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0),
                                       Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0),
                                       Kratos::make_intrusive<Node>(3, 0.0, 1.0, 0.0));

    const LinearNodalExtrapolator nodal_extrapolator;

    constexpr auto integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;
    auto           extrapolation_matrix =
        nodal_extrapolator.CalculateElementExtrapolationMatrix(geometry, integration_method);

    // clang-format off
    Matrix expected_extrapolation_matrix = ZeroMatrix(3, 3);
    expected_extrapolation_matrix <<= 5.0/3.0, -1.0/3.0, -1.0/3.0,
                                     -1.0/3.0,  5.0/3.0, -1.0/3.0,
                                     -1.0/3.0, -1.0/3.0,  5.0/3.0;
    // clang-format on

    KRATOS_EXPECT_MATRIX_NEAR(extrapolation_matrix, expected_extrapolation_matrix, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(NodalExtrapolator_GivesCorrectExtrapolationMatrix_For2D6NTriangle,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Kratos::Triangle2D6<Node> geometry(
        Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0), Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0),
        Kratos::make_intrusive<Node>(3, 0.0, 1.0, 0.0), Kratos::make_intrusive<Node>(4, 0.5, 0.0, 0.0),
        Kratos::make_intrusive<Node>(5, 0.5, 0.5, 0.0), Kratos::make_intrusive<Node>(6, 0.0, 0.5, 0.0));

    const LinearNodalExtrapolator nodal_extrapolator;

    constexpr auto integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;
    auto           extrapolation_matrix =
        nodal_extrapolator.CalculateElementExtrapolationMatrix(geometry, integration_method);

    // clang-format off
    Matrix expected_extrapolation_matrix = ZeroMatrix(6, 3);
    expected_extrapolation_matrix <<= 5.0/3.0, -1.0/3.0, -1.0/3.0,
                                     -1.0/3.0,  5.0/3.0, -1.0/3.0,
                                     -1.0/3.0, -1.0/3.0,  5.0/3.0,
                                      2.0/3.0,  2.0/3.0, -1.0/3.0,
                                     -1.0/3.0,  2.0/3.0,  2.0/3.0,
                                      2.0/3.0, -1.0/3.0,  2.0/3.0;
    // clang-format on

    KRATOS_EXPECT_MATRIX_NEAR(extrapolation_matrix, expected_extrapolation_matrix, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(NodalExtrapolator_GivesCorrectExtrapolationMatrix_For2D4NQuadrilateral,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Kratos::Quadrilateral2D4<Node> geometry(
        Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0), Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0),
        Kratos::make_intrusive<Node>(3, 1.0, 1.0, 0.0), Kratos::make_intrusive<Node>(4, 0.0, 1.0, 0.0));

    const LinearNodalExtrapolator nodal_extrapolator;

    constexpr auto integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;
    auto           extrapolation_matrix =
        nodal_extrapolator.CalculateElementExtrapolationMatrix(geometry, integration_method);

    // clang-format off
    Matrix expected_extrapolation_matrix = ZeroMatrix(4, 4);
    expected_extrapolation_matrix <<= 1.866025, -0.5,       0.133974, -0.5,
                                     -0.5,       1.866025, -0.5,       0.133974,
                                      0.133974, -0.5,       1.866025, -0.5,
                                     -0.5,       0.133974, -0.5,       1.866025;
    // clang-format on

    KRATOS_EXPECT_MATRIX_NEAR(extrapolation_matrix, expected_extrapolation_matrix, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(NodalExtrapolator_GivesCorrectExtrapolationMatrix_For2D8NQuadrilateral,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Kratos::Quadrilateral2D8<Node> geometry(
        Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0), Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0),
        Kratos::make_intrusive<Node>(3, 1.0, 1.0, 0.0), Kratos::make_intrusive<Node>(4, 0.0, 1.0, 0.0),
        Kratos::make_intrusive<Node>(5, 0.5, 0.0, 0.0), Kratos::make_intrusive<Node>(6, 1.0, 0.5, 0.0),
        Kratos::make_intrusive<Node>(7, 0.5, 1.0, 0.0), Kratos::make_intrusive<Node>(8, 0.0, 0.5, 0.0));

    const LinearNodalExtrapolator nodal_extrapolator;

    constexpr auto integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;
    auto           extrapolation_matrix =
        nodal_extrapolator.CalculateElementExtrapolationMatrix(geometry, integration_method);

    // clang-format off
    Matrix expected_extrapolation_matrix = ZeroMatrix(8, 4);
    expected_extrapolation_matrix <<= 1.866025, -0.5,       0.133974, -0.5,
                                     -0.5,       1.866025, -0.5,       0.133974,
                                      0.133974, -0.5,       1.866025, -0.5,
                                     -0.5,       0.133974, -0.5,       1.866025,
                                      0.683013,  0.683013, -0.183013, -0.183013,
                                     -0.183013,  0.683013,  0.683013, -0.183013,
                                     -0.183013, -0.183013,  0.683013,  0.683013,
                                      0.683013, -0.183013, -0.183013,  0.683013;
    // clang-format on

    KRATOS_EXPECT_MATRIX_NEAR(extrapolation_matrix, expected_extrapolation_matrix, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(NodalExtrapolator_GivesCorrectExtrapolationMatrix_For3D4NKratos_Tetrahedra,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Kratos::Tetrahedra3D4<Node> geometry(
        Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0), Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0),
        Kratos::make_intrusive<Node>(3, 0.0, 1.0, 0.0), Kratos::make_intrusive<Node>(4, 0.0, 0.0, 1.0));

    const LinearNodalExtrapolator nodal_extrapolator;

    constexpr auto integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;
    auto           extrapolation_matrix =
        nodal_extrapolator.CalculateElementExtrapolationMatrix(geometry, integration_method);

    // clang-format off
    Matrix expected_extrapolation_matrix = ZeroMatrix(4, 4);
    expected_extrapolation_matrix <<= -0.309017, -0.309017, -0.309017,  1.927051,
                                       1.927051, -0.309017, -0.309017, -0.309017,
                                      -0.309017,  1.927051, -0.309017, -0.309017,
                                      -0.309017, -0.309017,  1.927051, -0.309017;
    // clang-format on

    KRATOS_EXPECT_MATRIX_NEAR(extrapolation_matrix, expected_extrapolation_matrix, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(NodalExtrapolator_GivesCorrectExtrapolationMatrix_For3D10NKratos_Tetrahedra,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Kratos::Tetrahedra3D10<Node> geometry(
        Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0), Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0),
        Kratos::make_intrusive<Node>(3, 0.0, 1.0, 0.0), Kratos::make_intrusive<Node>(4, 0.0, 0.0, 1.0),
        Kratos::make_intrusive<Node>(5, 0.5, 0.0, 0.0), Kratos::make_intrusive<Node>(6, 0.5, 0.5, 0.0),
        Kratos::make_intrusive<Node>(7, 0.0, 0.5, 0.0), Kratos::make_intrusive<Node>(8, 0.0, 0.0, 0.5),
        Kratos::make_intrusive<Node>(9, 0.5, 0.0, 0.5), Kratos::make_intrusive<Node>(10, 0.0, 0.5, 0.5));

    const LinearNodalExtrapolator nodal_extrapolator;

    constexpr auto integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;
    auto           extrapolation_matrix =
        nodal_extrapolator.CalculateElementExtrapolationMatrix(geometry, integration_method);

    // clang-format off
    Matrix expected_extrapolation_matrix = ZeroMatrix(10, 4);
    expected_extrapolation_matrix <<= -0.309017, -0.309017, -0.309017,  1.927051,
                                       1.927051, -0.309017, -0.309017, -0.309017,
                                      -0.309017,  1.927051, -0.309017, -0.309017,
                                      -0.309017, -0.309017,  1.927051, -0.309017,

                                       0.809017, -0.309017, -0.309017,  0.809017,
                                       0.809017,  0.809017, -0.309017, -0.309017,
                                      -0.309017,  0.809017, -0.309017,  0.809017,
                                      -0.309017, -0.309017,  0.809017,  0.809017,
                                       0.809017, -0.309017,  0.809017, -0.309017,
                                      -0.309017,  0.809017,  0.809017, -0.309017;
    // clang-format on

    KRATOS_EXPECT_MATRIX_NEAR(extrapolation_matrix, expected_extrapolation_matrix, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(NodalExtrapolator_GivesCorrectExtrapolationMatrix_For3D8NKratos_Hexahedra,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Kratos::Hexahedra3D8<Node> geometry(
        Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0), Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0),
        Kratos::make_intrusive<Node>(3, 1.0, 1.0, 0.0), Kratos::make_intrusive<Node>(4, 0.0, 1.0, 0.0),
        Kratos::make_intrusive<Node>(5, 0.0, 0.0, 1.0), Kratos::make_intrusive<Node>(6, 1.0, 0.0, 1.0),
        Kratos::make_intrusive<Node>(7, 1.0, 1.0, 1.0), Kratos::make_intrusive<Node>(8, 0.0, 1.0, 1.0));

    const LinearNodalExtrapolator nodal_extrapolator;

    constexpr auto integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;
    auto           extrapolation_matrix =
        nodal_extrapolator.CalculateElementExtrapolationMatrix(geometry, integration_method);

    // clang-format off
    Matrix expected_extrapolation_matrix = ZeroMatrix(8, 8);
    expected_extrapolation_matrix <<=  2.549038, -0.683013,  0.183013, -0.683013, -0.683013,  0.183013, -0.049038,  0.183013,
                                      -0.683013,  2.549038, -0.683013,  0.183013,  0.183013, -0.683013,  0.183013, -0.049038,
                                       0.183013, -0.683013,  2.549038, -0.683013, -0.049038,  0.183013, -0.683013,  0.183013,
                                      -0.683013,  0.183013, -0.683013,  2.549038,  0.183013, -0.049038,  0.183013, -0.683013,
                                      -0.683013,  0.183013, -0.049038,  0.183013,  2.549038, -0.683013,  0.183013, -0.683013,
                                       0.183013, -0.683013,  0.183013, -0.049038, -0.683013,  2.549038, -0.683013,  0.183013,
                                      -0.049038,  0.183013, -0.683013,  0.183013,  0.183013, -0.683013,  2.549038, -0.683013,
                                       0.183013, -0.049038,  0.183013, -0.683013, -0.683013,  0.183013, -0.683013,  2.549038;
    // clang-format on

    KRATOS_EXPECT_MATRIX_NEAR(extrapolation_matrix, expected_extrapolation_matrix, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(NodalExtrapolator_GivesCorrectExtrapolationMatrix_For3D20NKratos_Hexahedra,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Kratos::Hexahedra3D20<Node> geometry(
        Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0), Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0),
        Kratos::make_intrusive<Node>(3, 1.0, 1.0, 0.0), Kratos::make_intrusive<Node>(4, 0.0, 1.0, 0.0),
        Kratos::make_intrusive<Node>(5, 0.0, 0.0, 1.0), Kratos::make_intrusive<Node>(6, 1.0, 0.0, 1.0),
        Kratos::make_intrusive<Node>(7, 1.0, 1.0, 1.0), Kratos::make_intrusive<Node>(8, 0.0, 1.0, 1.0),

        Kratos::make_intrusive<Node>(9, 0.5, 0.0, 0.0), Kratos::make_intrusive<Node>(10, 1.0, 0.5, 0.0),
        Kratos::make_intrusive<Node>(11, 0.5, 1.0, 0.0), Kratos::make_intrusive<Node>(12, 0.0, 0.5, 0.0),

        Kratos::make_intrusive<Node>(13, 0.0, 0.0, 0.5), Kratos::make_intrusive<Node>(14, 1.0, 0.0, 0.5),
        Kratos::make_intrusive<Node>(15, 1.0, 1.0, 0.5), Kratos::make_intrusive<Node>(16, 0.5, 1.0, 1.0),

        Kratos::make_intrusive<Node>(17, 0.5, 0.0, 1.0), Kratos::make_intrusive<Node>(18, 1.0, 0.5, 1.0),
        Kratos::make_intrusive<Node>(19, 0.5, 1.0, 1.0), Kratos::make_intrusive<Node>(20, 0.0, 0.5, 1.0));

    const LinearNodalExtrapolator nodal_extrapolator;

    constexpr auto integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;
    auto           extrapolation_matrix =
        nodal_extrapolator.CalculateElementExtrapolationMatrix(geometry, integration_method);

    // clang-format off
    Matrix expected_extrapolation_matrix = ZeroMatrix(20, 8);
    expected_extrapolation_matrix <<= 2.549038,  -0.683013,   0.183013,  -0.683013,  -0.683013,   0.183013,  -0.049038,   0.183013,
                                     -0.683013,   2.549038,  -0.683013,   0.183013,   0.183013,  -0.683013,   0.183013,  -0.049038,
                                      0.183013,  -0.683013,   2.549038,  -0.683013,  -0.049038,   0.183013,  -0.683013,   0.183013,
                                     -0.683013,   0.183013,  -0.683013,   2.549038,   0.183013,  -0.049038,   0.183013,  -0.683013,
                                     -0.683013,   0.183013,  -0.049038,   0.183013,   2.549038,  -0.683013,   0.183013,  -0.683013,
                                      0.183013,  -0.683013,   0.183013,  -0.049038,  -0.683013,   2.549038,  -0.683013,   0.183013,
                                     -0.049038,   0.183013,  -0.683013,   0.183013,   0.183013,  -0.683013,   2.549038,  -0.683013,
                                      0.183013,  -0.049038,   0.183013,  -0.683013,  -0.683013,   0.183013,  -0.683013,   2.549038,

                                      0.9330125,  0.9330125, -0.25,      -0.25,      -0.25,      -0.25,       0.0669875,  0.0669875,
                                     -0.25,       0.9330125,  0.9330125, -0.25,       0.0669875, -0.25,      -0.25,       0.0669875,
                                     -0.25,      -0.25,       0.9330125,  0.9330125,  0.0669875,  0.0669875, -0.25,      -0.25,
                                      0.9330125, -0.25,      -0.25,       0.9330125, -0.25,       0.0669875,  0.0669875, -0.25,
                                      0.9330125, -0.25,       0.0669875, -0.25,       0.9330125, -0.25,       0.0669875, -0.25,
                                     -0.25,       0.9330125, -0.25,       0.0669875, -0.25,       0.9330125, -0.25,       0.0669875,
                                      0.0669875, -0.25,       0.9330125, -0.25,       0.0669875, -0.25,       0.9330125, -0.25,
                                     -0.25,       0.0669875, -0.25,       0.9330125, -0.25,       0.0669875, -0.25,       0.9330125,
                                     -0.25,      -0.25,       0.0669875,  0.0669875,  0.9330125,  0.9330125, -0.25,      -0.25,
                                      0.0669875, -0.25,      -0.25,       0.0669875, -0.25,       0.9330125,  0.9330125, -0.25,
                                      0.0669875,  0.0669875, -0.25,      -0.25,      -0.25,      -0.25,       0.9330125,  0.9330125,
                                     -0.25,       0.0669875,  0.0669875, -0.25,       0.9330125, -0.25,      -0.25,       0.9330125;

    // clang-format on

    KRATOS_EXPECT_MATRIX_NEAR(extrapolation_matrix, expected_extrapolation_matrix, 1e-6)
}

} // namespace Kratos::Testing
