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
//                   Anne van de Graaf
//

#include "custom_geometries/line_interface_geometry.h"
#include "geometries/geometry_data.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_2d_3.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

#include <boost/numeric/ublas/assignment.hpp>

namespace
{

using namespace Kratos;

auto CreateTwoPlusTwoNoded2DLineInterfaceGeometry()
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 5.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, -1.0, 0.2, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 7.0, 0.2, 0.0));
    return LineInterfaceGeometry<Line2D2<Node>>{1, nodes};
}

auto CreateThreePlusThreeNoded2DLineInterfaceGeometry()
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 5.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, 2.5, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, -1.0, 0.2, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(5, 7.0, 0.2, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(6, 3.5, 0.4, 0.0));
    return LineInterfaceGeometry<Line2D3<Node>>{1, nodes};
}

} // namespace

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometryIsAGeometry, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry      = LineInterfaceGeometry<Line2D3<Node>>();
    const auto base_geometry = dynamic_cast<const Geometry<Node>*>(&geometry);

    KRATOS_EXPECT_NE(base_geometry, nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometryCanBeConstructedGivenASetOfNullPointersToNodes,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // The following constructor input data resembles what is done at element registration time
    const auto six_null_pointers_to_nodes = Geometry<Node>::PointsArrayType{6};

    const auto geometry = LineInterfaceGeometry<Line2D3<Node>>{six_null_pointers_to_nodes};

    KRATOS_EXPECT_EQ(geometry.PointsNumber(), 6);
    KRATOS_EXPECT_EQ(geometry.LocalSpaceDimension(), 1);
    KRATOS_EXPECT_EQ(geometry.WorkingSpaceDimension(), 2);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_Create_CreatesNewInstanceOfCorrectType, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto          geometry = LineInterfaceGeometry<Line2D2<Node>>();
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 0.0, 0.0, 0.0));

    const auto new_geometry = geometry.Create(nodes);

    KRATOS_EXPECT_NE(new_geometry, nullptr);
    KRATOS_EXPECT_NE(dynamic_cast<const LineInterfaceGeometry<Line2D2<Node>>*>(new_geometry.get()), nullptr);
    KRATOS_EXPECT_EQ(new_geometry->PointsNumber(), 4);
    KRATOS_EXPECT_EQ(new_geometry->Id(), 0);
    KRATOS_EXPECT_EQ(new_geometry->LocalSpaceDimension(), 1);
    KRATOS_EXPECT_EQ(new_geometry->WorkingSpaceDimension(), 2);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_CreateWithId_CreatesNewInstanceOfCorrectTypeAndId,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto          geometry = LineInterfaceGeometry<Line2D2<Node>>();
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 0.0, 0.0, 0.0));

    constexpr auto new_geometry_id = 1;
    const auto     new_geometry    = geometry.Create(new_geometry_id, nodes);

    KRATOS_EXPECT_NE(new_geometry, nullptr);
    KRATOS_EXPECT_NE(dynamic_cast<const LineInterfaceGeometry<Line2D2<Node>>*>(new_geometry.get()), nullptr);
    KRATOS_EXPECT_EQ(new_geometry->PointsNumber(), 4);
    KRATOS_EXPECT_EQ(new_geometry->Id(), new_geometry_id);
    KRATOS_EXPECT_EQ(new_geometry->LocalSpaceDimension(), 1);
    KRATOS_EXPECT_EQ(new_geometry->WorkingSpaceDimension(), 2);
}

KRATOS_TEST_CASE_IN_SUITE(CreatingInterfaceWithThreeNodesThrows, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 5.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, -1.0, 0.0, 0.0));

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(LineInterfaceGeometry<Line2D3<Node>>{nodes},
                                      "Number of nodes must be 2+2 or 3+3")
    constexpr auto geometry_id = 1;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN((LineInterfaceGeometry<Line2D3<Node>>{geometry_id, nodes}),
                                      "Number of nodes must be 2+2 or 3+3")
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectShapeFunctionValuesInNodes_ForTwoPlusTwoNodedGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateTwoPlusTwoNoded2DLineInterfaceGeometry();
    const auto xi_start = array_1d<double, 3>{-1.0, 0.0, 0.0};
    const auto xi_end   = array_1d<double, 3>{1.0, 0.0, 0.0};

    // Values for first shape function
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(0, xi_start), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(0, xi_end), 0.0);

    // Values for second shape function
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(1, xi_start), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(1, xi_end), 1.0);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectShapeFunctionValuesInNodes_ForThreePlusThreeNodedGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry  = CreateThreePlusThreeNoded2DLineInterfaceGeometry();
    const auto xi_start  = array_1d<double, 3>{-1.0, 0.0, 0.0};
    const auto xi_end    = array_1d<double, 3>{1.0, 0.0, 0.0};
    const auto xi_middle = array_1d<double, 3>{0.0, 0.0, 0.0};

    // Values for first shape function
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(0, xi_start), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(0, xi_end), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(0, xi_middle), 0.0);

    // Values for second shape function
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(1, xi_start), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(1, xi_end), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(1, xi_middle), 0.0);

    // Values for third shape function
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(2, xi_start), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(2, xi_end), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(2, xi_middle), 1.0);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectAllShapeFunctionValuesAtPosition_ForTwoPlusTwoNodedGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateTwoPlusTwoNoded2DLineInterfaceGeometry();
    const auto xi       = array_1d<double, 3>{0.5, 0.0, 0.0};

    Vector result;
    geometry.ShapeFunctionsValues(result, xi);

    // Note that the shape function values are evaluated per nodal pair!
    Vector expected_result{2};
    expected_result <<= 0.25, 0.75;
    KRATOS_EXPECT_VECTOR_NEAR(result, expected_result, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectAllShapeFunctionValuesAtPosition_ForThreePlusThreeNodedGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateThreePlusThreeNoded2DLineInterfaceGeometry();
    const auto xi       = array_1d<double, 3>{0.5, 0.0, 0.0};

    Vector result;
    geometry.ShapeFunctionsValues(result, xi);

    // Note that the shape function values are evaluated per nodal pair!
    Vector expected_result{3};
    expected_result <<= -0.125, 0.375, 0.75;
    KRATOS_EXPECT_VECTOR_NEAR(result, expected_result, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectAllLocalGradientsAtPosition_ForTwoPlusTwoNodedGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateTwoPlusTwoNoded2DLineInterfaceGeometry();
    const auto xi       = array_1d<double, 3>{0.5, 0.0, 0.0};

    Matrix result;
    geometry.ShapeFunctionsLocalGradients(result, xi);

    Matrix expected_result(2, 1);
    expected_result <<= -0.5, 0.5;
    KRATOS_EXPECT_MATRIX_NEAR(result, expected_result, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectAllLocalGradientsAtPosition_ForThreePlusThreeNodedGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateThreePlusThreeNoded2DLineInterfaceGeometry();
    const auto xi       = array_1d<double, 3>{0.5, 0.0, 0.0};

    Matrix result;
    geometry.ShapeFunctionsLocalGradients(result, xi);

    Matrix expected_result(3, 1);
    expected_result <<= 0.0, 1.0, -1.0;
    KRATOS_EXPECT_MATRIX_NEAR(result, expected_result, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectJacobian_ForTwoPlusTwoNodedGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateTwoPlusTwoNoded2DLineInterfaceGeometry();
    const auto xi       = array_1d<double, 3>{0.5, 0.0, 0.0};

    Matrix result;
    geometry.Jacobian(result, xi);

    Matrix expected_result(2, 1);
    expected_result <<= 3.25, 0.0;
    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(result, expected_result, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectJacobian_ForThreePlusThreeNodedGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateThreePlusThreeNoded2DLineInterfaceGeometry();
    const auto xi       = array_1d<double, 3>{0.5, 0.0, 0.0};

    Matrix result;
    geometry.Jacobian(result, xi);

    Matrix expected_result(2, 1);
    expected_result <<= 3.0, -0.1;
    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(result, expected_result, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectDeterminantOfJacobian_ForTwoPlusTwoNodedGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateTwoPlusTwoNoded2DLineInterfaceGeometry();
    const auto xi       = array_1d<double, 3>{0.5, 0.0, 0.0};

    KRATOS_EXPECT_RELATIVE_NEAR(geometry.DeterminantOfJacobian(xi), 3.25, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectDeterminantOfJacobian_ForThreePlusThreeNodedGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateThreePlusThreeNoded2DLineInterfaceGeometry();
    const auto xi       = array_1d<double, 3>{0.5, 0.0, 0.0};

    KRATOS_EXPECT_RELATIVE_NEAR(geometry.DeterminantOfJacobian(xi), std::sqrt(9.01), 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectWorkingSpaceDimension, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateThreePlusThreeNoded2DLineInterfaceGeometry();
    KRATOS_EXPECT_EQ(geometry.WorkingSpaceDimension(), 2);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_Throws_WhenCallingInverseJacobian, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateTwoPlusTwoNoded2DLineInterfaceGeometry();
    const auto xi       = array_1d<double, 3>{0.5, 0.0, 0.0};

    Matrix result;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        geometry.InverseOfJacobian(result, xi),
        "Inverse of Jacobian is not implemented for the line interface geometry")
}

KRATOS_TEST_CASE_IN_SUITE(TwoPlusTwoLineInterfaceGeometry_LengthReturnsTheLengthOfUnderlyingLineGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateTwoPlusTwoNoded2DLineInterfaceGeometry();

    constexpr auto expected_length = 6.5;
    KRATOS_EXPECT_RELATIVE_NEAR(geometry.Length(), expected_length, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(ThreePlusThreeNodedLineInterfaceGeometry_LengthReturnsTheLengthOfUnderlyingLineGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateThreePlusThreeNoded2DLineInterfaceGeometry();

    // This number was not calculated by hand, meaning this unit test is a regression test.
    constexpr auto expected_length = 6.504159;
    KRATOS_EXPECT_RELATIVE_NEAR(geometry.Length(), expected_length, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(TwoPlusTwoNodedLineInterfaceGeometry_DomainSizeReturnsTheLengthOfUnderlyingLineGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateTwoPlusTwoNoded2DLineInterfaceGeometry();

    KRATOS_EXPECT_RELATIVE_NEAR(geometry.DomainSize(), geometry.Length(), 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(ThreePlusThreeNodedLineInterfaceGeometry_DomainSizeReturnsTheLengthOfUnderlyingLineGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateThreePlusThreeNoded2DLineInterfaceGeometry();

    KRATOS_EXPECT_RELATIVE_NEAR(geometry.DomainSize(), geometry.Length(), 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(GlobalCoordinatesAreCorrectlyMappedToLocalCoordinate_ForTwoPlusTwoNodedLineInterfaceGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateTwoPlusTwoNoded2DLineInterfaceGeometry();

    const auto global_position = array_1d<double, 3>{4.375, 0.1, 0.0};
    auto       result          = array_1d<double, 3>{};
    result                     = geometry.PointLocalCoordinates(result, global_position);

    const auto expected_result = array_1d<double, 3>{0.5, 0.0, 0.0};
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(result, expected_result, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(GlobalCoordinatesAreCorrectlyMappedToLocalCoordinate_ForThreePlusThreeNodedLineInterfaceGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateThreePlusThreeNoded2DLineInterfaceGeometry();

    const auto global_position = array_1d<double, 3>{4.5625, 0.175, 0.0};
    auto       result          = array_1d<double, 3>{};
    result                     = geometry.PointLocalCoordinates(result, global_position);

    const auto expected_result = array_1d<double, 3>{0.5, 0.0, 0.0};
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(result, expected_result, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(GetLocalCoordinatesOfAllNodesOfThreePlusThreeNodedLineInterfaceGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateThreePlusThreeNoded2DLineInterfaceGeometry();

    Matrix result;
    geometry.PointsLocalCoordinates(result);

    Matrix expected_result{3, 1};
    expected_result <<= -1.0, 1.0, 0.0;
    KRATOS_EXPECT_MATRIX_NEAR(result, expected_result, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_Throws_WhenCallingArea, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateThreePlusThreeNoded2DLineInterfaceGeometry();
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(geometry.Area(),
                                      "Calling base class 'Area' method instead of derived class "
                                      "one. Please check the definition of derived class. ")
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_Throws_WhenCallingVolume, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateThreePlusThreeNoded2DLineInterfaceGeometry();
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(geometry.Volume(),
                                      "Calling base class 'Volume' method instead of derived class "
                                      "one. Please check the definition of derived class. ")
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_Throws_WhenCallingFunctionsRelatedToIntegrationPoints,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto                          geometry = CreateThreePlusThreeNoded2DLineInterfaceGeometry();
    Geometry<Node>::JacobiansType dummy_jacobian;
    Geometry<Node>::ShapeFunctionsGradientsType dummy_shape_functions_gradients;
    Vector                                      dummy_vector;
    Matrix                                      dummy_matrix;
    IndexType                                   dummy_index = 0;
    GeometryData::IntegrationMethod dummy_integration_method = GeometryData::IntegrationMethod::GI_GAUSS_1;
    Geometry<Node>::IntegrationPointsArrayType dummy_integration_points;
    IntegrationInfo dummy_integration_info(dummy_index, dummy_integration_method);
    Geometry<Node>::GeometriesArrayType     dummy_geometries;
    std::vector<Node::CoordinatesArrayType> dummy_coordinates;
    const std::string                       message =
        "This Geometry type does not support functionality related to integration schemes.\n";

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(geometry.Normal(dummy_index), message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(geometry.Normal(dummy_index, dummy_integration_method), message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(geometry.UnitNormal(dummy_index), message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(geometry.UnitNormal(dummy_index, dummy_integration_method), message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(geometry.Jacobian(dummy_jacobian, dummy_integration_method), message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        geometry.Jacobian(dummy_jacobian, dummy_integration_method, dummy_matrix), message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        geometry.Jacobian(dummy_matrix, dummy_index, dummy_integration_method), message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        geometry.Jacobian(dummy_matrix, dummy_index, dummy_integration_method, dummy_matrix), message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        geometry.DeterminantOfJacobian(dummy_vector, dummy_integration_method), message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(geometry.DeterminantOfJacobian(dummy_index, dummy_integration_method), message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(geometry.InverseOfJacobian(dummy_jacobian, dummy_integration_method), message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        geometry.InverseOfJacobian(dummy_matrix, dummy_index, dummy_integration_method), message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(geometry.ShapeFunctionsIntegrationPointsGradients(
                                          dummy_shape_functions_gradients, dummy_integration_method),
                                      message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(geometry.ShapeFunctionsIntegrationPointsGradients(
                                          dummy_shape_functions_gradients, dummy_vector, dummy_integration_method),
                                      message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        geometry.ShapeFunctionsIntegrationPointsGradients(
            dummy_shape_functions_gradients, dummy_vector, dummy_integration_method, dummy_matrix),
        message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        geometry.CreateIntegrationPoints(dummy_integration_points, dummy_integration_info), message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(geometry.GetDefaultIntegrationInfo(), message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        geometry.CreateQuadraturePointGeometries(dummy_geometries, dummy_index,
                                                 dummy_integration_points, dummy_integration_info),
        message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        geometry.CreateQuadraturePointGeometries(dummy_geometries, dummy_index, dummy_integration_info), message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        geometry.GlobalSpaceDerivatives(dummy_coordinates, dummy_index, dummy_index), message)
}
} // namespace Kratos::Testing