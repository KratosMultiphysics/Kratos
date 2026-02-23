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
#include "tests/cpp_tests/test_utilities.h"

using namespace Kratos;

namespace
{
struct Line2D2NGeometryWithVelocityVectorData {
    Model                          model;
    std::unique_ptr<Line2D2<Node>> p_geometry;
    array_1d<double, 3>            velocity_vector;

    Line2D2NGeometryWithVelocityVectorData()
    {
        ModelPart& r_model_part = model.CreateModelPart("Main");
        r_model_part.AddNodalSolutionStepVariable(VELOCITY);
        PointerVector<Node> nodes;
        nodes.push_back(r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0));
        nodes.push_back(r_model_part.CreateNewNode(2, 5.0, 0.0, 0.0));
        p_geometry      = std::make_unique<Line2D2<Node>>(nodes);
        velocity_vector = array_1d<double, 3>{1.0, 0.5, -0.5};
        for (auto& node : *p_geometry) {
            node.FastGetSolutionStepValue(VELOCITY) = velocity_vector;
        }
    }

    Vector ExpectedVelocity2DFlat() const
    {
        auto result            = Vector(4);
        subrange(result, 0, 2) = subrange(velocity_vector, 0, 2);
        subrange(result, 2, 4) = subrange(velocity_vector, 0, 2);
        return result;
    }

    Vector ExpectedVelocity3DFlat() const
    {
        auto result            = Vector(6);
        subrange(result, 0, 3) = velocity_vector;
        subrange(result, 3, 6) = velocity_vector;
        return result;
    }
};
} // namespace

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

KRATOS_TEST_CASE_IN_SUITE(ElementUtilities_GetNodalVariableVector_Returns2DFlatVector, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto     test_data       = Line2D2NGeometryWithVelocityVectorData();
    constexpr auto dimension_2D    = std::size_t{2};
    const auto     number_of_nodes = test_data.p_geometry->PointsNumber();

    // Act
    const auto actual_values_in_2D = GeoElementUtilities::GetNodalVariableVector(
        *test_data.p_geometry, VELOCITY, dimension_2D, dimension_2D * number_of_nodes);

    // Assert
    KRATOS_EXPECT_VECTOR_NEAR(actual_values_in_2D, test_data.ExpectedVelocity2DFlat(),
                              Testing::Defaults::relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(ElementUtilities_GetNodalVariableVector_Returns3DFlatVector, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto     test_data       = Line2D2NGeometryWithVelocityVectorData();
    constexpr auto dimension_3D    = std::size_t{3};
    const auto     number_of_nodes = test_data.p_geometry->PointsNumber();

    // Act
    const auto actual_values_in_2D = GeoElementUtilities::GetNodalVariableVector(
        *test_data.p_geometry, VELOCITY, dimension_3D, dimension_3D * number_of_nodes);

    // Assert
    KRATOS_EXPECT_VECTOR_NEAR(actual_values_in_2D, test_data.ExpectedVelocity3DFlat(),
                              Testing::Defaults::relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(ElementUtilities_GetNodalVariableVector_Returns2DFlatArray, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto     test_data       = Line2D2NGeometryWithVelocityVectorData();
    constexpr auto dimension_2D    = std::size_t{2};
    const auto     number_of_nodes = std::size_t{2};

    // Act
    auto actual_values_in_2D = array_1d<double, dimension_2D * number_of_nodes>();
    GeoElementUtilities::GetNodalVariableVector<dimension_2D, number_of_nodes>(
        actual_values_in_2D, *test_data.p_geometry, VELOCITY);

    // Assert
    ASSERT_EQ(number_of_nodes, test_data.p_geometry->PointsNumber());
    KRATOS_EXPECT_VECTOR_NEAR(actual_values_in_2D, test_data.ExpectedVelocity2DFlat(),
                              Testing::Defaults::relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(ElementUtilities_GetNodalVariableVector_Returns3DFlatArray, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto     test_data       = Line2D2NGeometryWithVelocityVectorData();
    constexpr auto dimension_3D    = std::size_t{3};
    const auto     number_of_nodes = std::size_t{2};

    // Act
    auto actual_values_in_3D = array_1d<double, dimension_3D * number_of_nodes>();
    GeoElementUtilities::GetNodalVariableVector<dimension_3D, number_of_nodes>(
        actual_values_in_3D, *test_data.p_geometry, VELOCITY);

    // Assert
    ASSERT_EQ(number_of_nodes, test_data.p_geometry->PointsNumber());
    KRATOS_EXPECT_VECTOR_NEAR(actual_values_in_3D, test_data.ExpectedVelocity3DFlat(),
                              Testing::Defaults::relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(ElementUtilities_GetNodalVariableVector_ChecksDimension, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto     test_data       = Line2D2NGeometryWithVelocityVectorData();
    constexpr auto wrong_dimension = std::size_t{4};
    const auto     number_of_nodes = test_data.p_geometry->PointsNumber();
    // Act and Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(GeoElementUtilities::GetNodalVariableVector(
                                          *test_data.p_geometry, VELOCITY, wrong_dimension, number_of_nodes),
                                      "Incorrect dimension value (4).");
}

KRATOS_TEST_CASE_IN_SUITE(ElementUtilities_GetNodalVariableVector_ChecksDOFs, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto     test_data       = Line2D2NGeometryWithVelocityVectorData();
    constexpr auto dimension_2D    = std::size_t{2};
    const auto     number_of_nodes = test_data.p_geometry->PointsNumber();
    // Act and Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        GeoElementUtilities::GetNodalVariableVector(*test_data.p_geometry, VELOCITY, dimension_2D, number_of_nodes), "Mismatch between requested number of DOFs (2) and computed number of DOFs from nodal values (4).");
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

class ParametrizedMatrixAssignAndAssembleFixture
    : public ::testing::TestWithParam<std::tuple<std::function<void(Matrix&, const Matrix&)>, Matrix>>
{
};

TEST_P(ParametrizedMatrixAssignAndAssembleFixture, MatricesAreFilledForVaryingParts)
{
    // Arrange
    const auto& [matrix_function, expected_matrix] = GetParam();

    auto       destination_matrix = Matrix{3, 3, 5.0};
    const auto addition_matrix    = Matrix{2, 2, 3.0};

    matrix_function(destination_matrix, addition_matrix);
    KRATOS_EXPECT_MATRIX_EQ(destination_matrix, expected_matrix);
}

INSTANTIATE_TEST_SUITE_P(
    KratosGeoMechanicsFastSuiteWithoutKernel,
    ParametrizedMatrixAssignAndAssembleFixture,
    ::testing::Values(
        std::make_tuple(GeoElementUtilities::AssignUUBlockMatrix<Matrix, Matrix>,
                        UblasUtilities::CreateMatrix({{3.0, 3.0, 5.0}, {3.0, 3.0, 5.0}, {5.0, 5.0, 5.0}})),
        std::make_tuple(GeoElementUtilities::AssignUPBlockMatrix<Matrix, Matrix>,
                        UblasUtilities::CreateMatrix({{5.0, 3.0, 3.0}, {5.0, 3.0, 3.0}, {5.0, 5.0, 5.0}})),
        std::make_tuple(GeoElementUtilities::AssignPUBlockMatrix<Matrix, Matrix>,
                        UblasUtilities::CreateMatrix({{5.0, 5.0, 5.0}, {3.0, 3.0, 5.0}, {3.0, 3.0, 5.0}})),
        std::make_tuple(GeoElementUtilities::AssignPPBlockMatrix<Matrix, Matrix>,
                        UblasUtilities::CreateMatrix({{5.0, 5.0, 5.0}, {5.0, 3.0, 3.0}, {5.0, 3.0, 3.0}})),
        std::make_tuple(GeoElementUtilities::AssembleUUBlockMatrix<Matrix, Matrix>,
                        UblasUtilities::CreateMatrix({{8.0, 8.0, 5.0}, {8.0, 8.0, 5.0}, {5.0, 5.0, 5.0}})),
        std::make_tuple(GeoElementUtilities::AssembleUPBlockMatrix<Matrix, Matrix>,
                        UblasUtilities::CreateMatrix({{5.0, 8.0, 8.0}, {5.0, 8.0, 8.0}, {5.0, 5.0, 5.0}})),
        std::make_tuple(GeoElementUtilities::AssemblePUBlockMatrix<Matrix, Matrix>,
                        UblasUtilities::CreateMatrix({{5.0, 5.0, 5.0}, {8.0, 8.0, 5.0}, {8.0, 8.0, 5.0}})),
        std::make_tuple(GeoElementUtilities::AssemblePPBlockMatrix<Matrix, Matrix>,
                        UblasUtilities::CreateMatrix({{5.0, 5.0, 5.0}, {5.0, 8.0, 8.0}, {5.0, 8.0, 8.0}}))));

} // namespace Kratos::Testing
