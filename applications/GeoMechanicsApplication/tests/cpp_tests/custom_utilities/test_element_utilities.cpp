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
//

#include "custom_geometries/interface_geometry.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/ublas_utilities.h"
#include "geometries/line_2d_2.h"
#include "includes/node.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(ElementUtilities_ReturnsCorrectListOfShapeFunctionsValuesAtIntegrationPoints,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 5.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, -1.0, 0.2, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 7.0, 0.2, 0.0));
    const auto geometry = InterfaceGeometry<Line2D2<Node>>{1, nodes};

    const Geo::IntegrationPointVectorType integration_points{{-1.0, 0.0, 0.0, 1.0}, {1.0, 0.0, 0.0, 1.0}};
    const auto shape_function_values =
        GeoElementUtilities::EvaluateShapeFunctionsAtIntegrationPoints(integration_points, geometry);

    auto shape_function_values_1 = UblasUtilities::CreateVector({1.0, 0.0});
    auto shape_function_values_2 = UblasUtilities::CreateVector({0.0, 1.0});
    const std::vector<Vector> expected_shape_function_values{shape_function_values_1, shape_function_values_2};

    KRATOS_EXPECT_EQ(expected_shape_function_values.size(), shape_function_values.size());
    for (std::size_t i = 0; i < expected_shape_function_values.size(); ++i) {
        KRATOS_CHECK_VECTOR_NEAR(expected_shape_function_values[i], shape_function_values[i], 1e-6);
    }
}

KRATOS_TEST_CASE_IN_SUITE(AssignUUBlockMatrix_PositionsMatrixAtTopLeft, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto destination_matrix = Matrix{3, 3, 5.0};
    auto insertion_matrix   = Matrix{2, 2, 3.0};

    // Act
    GeoElementUtilities::AssignUUBlockMatrix(destination_matrix, insertion_matrix);

    // Assert
    auto expected_matrix = UblasUtilities::CreateMatrix({{3.0, 3.0, 5.0}, {3.0, 3.0, 5.0}, {5.0, 5.0, 5.0}});
    KRATOS_EXPECT_MATRIX_EQ(destination_matrix, expected_matrix);
}

KRATOS_TEST_CASE_IN_SUITE(AssignUPBlockMatrix_PositionsMatrixAtTopRight, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto destination_matrix = Matrix{3, 3, 5.0};
    auto insertion_matrix   = Matrix{2, 2, 3.0};

    // Act
    GeoElementUtilities::AssignUPBlockMatrix(destination_matrix, insertion_matrix);

    // Assert
    auto expected_matrix = UblasUtilities::CreateMatrix({{5.0, 3.0, 3.0}, {5.0, 3.0, 3.0}, {5.0, 5.0, 5.0}});
    KRATOS_EXPECT_MATRIX_EQ(destination_matrix, expected_matrix);
}

KRATOS_TEST_CASE_IN_SUITE(AssignPUBlockMatrix_PositionsMatrixAtBottomLeft, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto destination_matrix = Matrix{3, 3, 5.0};
    auto insertion_matrix   = Matrix{2, 2, 3.0};

    // Act
    GeoElementUtilities::AssignPUBlockMatrix(destination_matrix, insertion_matrix);

    // Assert
    auto expected_matrix = UblasUtilities::CreateMatrix({{5.0, 5.0, 5.0}, {3.0, 3.0, 5.0}, {3.0, 3.0, 5.0}});
    KRATOS_EXPECT_MATRIX_EQ(destination_matrix, expected_matrix);
}

KRATOS_TEST_CASE_IN_SUITE(AssignPPBlockMatrix_PositionsMatrixAtBottomRight, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto destination_matrix = Matrix{3, 3, 5.0};
    auto insertion_matrix   = Matrix{2, 2, 3.0};

    // Act
    GeoElementUtilities::AssignPPBlockMatrix(destination_matrix, insertion_matrix);

    // Assert
    auto expected_matrix = UblasUtilities::CreateMatrix({{5.0, 5.0, 5.0}, {5.0, 3.0, 3.0}, {5.0, 3.0, 3.0}});
    KRATOS_EXPECT_MATRIX_EQ(destination_matrix, expected_matrix);
}

KRATOS_TEST_CASE_IN_SUITE(AssembleUUBlockMatrix_AddsMatrixAtTopLeft, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto destination_matrix = Matrix{3, 3, 5.0};
    auto addition_matrix    = Matrix{2, 2, 3.0};

    // Act
    GeoElementUtilities::AssembleUUBlockMatrix(destination_matrix, addition_matrix);

    // Assert
    auto expected_matrix = UblasUtilities::CreateMatrix({{8.0, 8.0, 5.0}, {8.0, 8.0, 5.0}, {5.0, 5.0, 5.0}});
    KRATOS_EXPECT_MATRIX_EQ(destination_matrix, expected_matrix);
}

KRATOS_TEST_CASE_IN_SUITE(AssembleUPBlockMatrix_AddsMatrixAtTopRight, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto destination_matrix = Matrix{3, 3, 5.0};
    auto addition_matrix    = Matrix{2, 2, 3.0};

    // Act
    GeoElementUtilities::AssembleUPBlockMatrix(destination_matrix, addition_matrix);

    // Assert
    auto expected_matrix = UblasUtilities::CreateMatrix({{5.0, 8.0, 8.0}, {5.0, 8.0, 8.0}, {5.0, 5.0, 5.0}});
    KRATOS_EXPECT_MATRIX_EQ(destination_matrix, expected_matrix);
}

KRATOS_TEST_CASE_IN_SUITE(AssemblePUBlockMatrix_AddsMatrixAtBottomLeft, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto destination_matrix = Matrix{3, 3, 5.0};
    auto addition_matrix    = Matrix{2, 2, 3.0};

    // Act
    GeoElementUtilities::AssemblePUBlockMatrix(destination_matrix, addition_matrix);

    // Assert
    auto expected_matrix = UblasUtilities::CreateMatrix({{5.0, 5.0, 5.0}, {8.0, 8.0, 5.0}, {8.0, 8.0, 5.0}});
    KRATOS_EXPECT_MATRIX_EQ(destination_matrix, expected_matrix);
}

KRATOS_TEST_CASE_IN_SUITE(AssemblePPBlockMatrix_AddsMatrixAtBottomRight, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto destination_matrix = Matrix{3, 3, 5.0};
    auto addition_matrix    = Matrix{2, 2, 3.0};

    // Act
    GeoElementUtilities::AssemblePPBlockMatrix(destination_matrix, addition_matrix);

    // Assert
    auto expected_matrix = UblasUtilities::CreateMatrix({{5.0, 5.0, 5.0}, {5.0, 8.0, 8.0}, {5.0, 8.0, 8.0}});
    KRATOS_EXPECT_MATRIX_EQ(destination_matrix, expected_matrix);
}

class ParametrizedVectorAssignAndAssembleFixture
    : public ::testing::TestWithParam<std::tuple<std::function<void(Vector&, const Vector&)>, Vector>>
{
};

TEST_P(ParametrizedVectorAssignAndAssembleFixture, VectorsAreFilledForVaryingParts)
{
    // Arrange
    const auto& [vector_function, expected_vector] = GetParam();

    auto       destination_vector = Vector{3, 5.0};
    const auto addition_vector    = Vector{2, 3.0};

    vector_function(destination_vector, addition_vector);
    KRATOS_EXPECT_VECTOR_EQ(destination_vector, expected_vector);
}

INSTANTIATE_TEST_SUITE_P(
    KratosGeoMechanicsFastSuiteWithoutKernel,
    ParametrizedVectorAssignAndAssembleFixture,
    ::testing::Values(std::make_tuple(GeoElementUtilities::AssignUBlockVector<Vector, Vector>,
                                      UblasUtilities::CreateVector({3.0, 3.0, 5.0})),
                      std::make_tuple(GeoElementUtilities::AssignPBlockVector<Vector, Vector>,
                                      UblasUtilities::CreateVector({5.0, 3.0, 3.0})),
                      std::make_tuple(GeoElementUtilities::AssembleUBlockVector<Vector, Vector>,
                                      UblasUtilities::CreateVector({8.0, 8.0, 5.0})),
                      std::make_tuple(GeoElementUtilities::AssemblePBlockVector<Vector, Vector>,
                                      UblasUtilities::CreateVector({5.0, 8.0, 8.0}))));
} // namespace Kratos::Testing
