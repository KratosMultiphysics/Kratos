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

#include "../geo_mechanics_fast_suite.h"
#include "custom_geometries/line_interface_geometry.h"
#include <boost/numeric/ublas/assignment.hpp>

namespace
{

using namespace Kratos;

LineInterfaceGeometry CreateFourNodedLineInterfaceGeometry()
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 5.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, -1.0, 0.2, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 7.0, 0.2, 0.0));
    return {1, nodes};
}

LineInterfaceGeometry CreateSixNodedLineInterfaceGeometry()
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 5.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, 2.5, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, -1.0, 0.2, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(5, 7.0, 0.2, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(6, 3.5, 0.4, 0.0));
    return {1, nodes};
}

} // namespace

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometryIsAGeometry, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry      = LineInterfaceGeometry();
    auto       base_geometry = dynamic_cast<const Geometry<Node>*>(&geometry);

    KRATOS_EXPECT_NE(base_geometry, nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_Create_CreatesNewInstanceOfCorrectType, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto          geometry = LineInterfaceGeometry();
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 0.0, 0.0, 0.0));

    const auto new_geometry = geometry.Create(nodes);

    KRATOS_EXPECT_NE(new_geometry, nullptr);
    KRATOS_EXPECT_NE(dynamic_cast<const LineInterfaceGeometry*>(new_geometry.get()), nullptr);
    KRATOS_EXPECT_EQ(new_geometry->PointsNumber(), 4);
    KRATOS_EXPECT_EQ(new_geometry->Id(), 0);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_CreateWithId_CreatesNewInstanceOfCorrectTypeAndId,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto          geometry = LineInterfaceGeometry();
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 0.0, 0.0, 0.0));

    const auto new_geometry = geometry.Create(1, nodes);

    KRATOS_EXPECT_NE(new_geometry, nullptr);
    KRATOS_EXPECT_NE(dynamic_cast<const LineInterfaceGeometry*>(new_geometry.get()), nullptr);
    KRATOS_EXPECT_EQ(new_geometry->PointsNumber(), 4);
    KRATOS_EXPECT_EQ(new_geometry->Id(), 1);
}

KRATOS_TEST_CASE_IN_SUITE(ExpectThrowWhenCreatingInterfaceGivenTwoNodes, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 0.0, 0.0, 0.0));

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(LineInterfaceGeometry{nodes},
                                      "Number of nodes must be four or six")
}

KRATOS_TEST_CASE_IN_SUITE(CreatingInterfaceWithOddNumberOfNodesThrows, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 5.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, -1.0, 0.0, 0.0));

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(LineInterfaceGeometry{nodes},
                                      "Number of nodes must be four or six")
    constexpr auto geometry_id = 1;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN((LineInterfaceGeometry{geometry_id, nodes}),
                                      "Number of nodes must be four or six")
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectShapeFunctionValuesInNodes_ForFourNodedGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateFourNodedLineInterfaceGeometry();

    const auto ksi_start = array_1d<double, 3>{-1.0, 0.0, 0.0};
    const auto ksi_end   = array_1d<double, 3>{1.0, 0.0, 0.0};

    // Values for first shape function
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(0, ksi_start), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(0, ksi_end), 0.0);

    // Values for second shape function
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(1, ksi_start), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(1, ksi_end), 1.0);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectShapeFunctionValuesInNodes_ForSixNodedGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateSixNodedLineInterfaceGeometry();

    const auto ksi_start  = array_1d<double, 3>{-1.0, 0.0, 0.0};
    const auto ksi_end    = array_1d<double, 3>{1.0, 0.0, 0.0};
    const auto ksi_middle = array_1d<double, 3>{0.0, 0.0, 0.0};

    // Values for first shape function
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(0, ksi_start), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(0, ksi_end), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(0, ksi_middle), 0.0);

    // Values for second shape function
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(1, ksi_start), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(1, ksi_end), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(1, ksi_middle), 0.0);

    // Values for third shape function
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(2, ksi_start), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(2, ksi_end), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(2, ksi_middle), 1.0);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectAllShapeFunctionValuesAtPosition_ForFourNodedGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateFourNodedLineInterfaceGeometry();

    const auto ksi = array_1d<double, 3>{0.5, 0.0, 0.0};

    Vector result;
    geometry.ShapeFunctionsValues(result, ksi);

    // Note that the shape function values are evaluated per nodal pair!
    Vector expected_result{2};
    expected_result <<= 0.25, 0.75;
    KRATOS_EXPECT_VECTOR_NEAR(result, expected_result, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectAllShapeFunctionValuesAtPosition_ForSixNodedGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateSixNodedLineInterfaceGeometry();

    const auto ksi = array_1d<double, 3>{0.5, 0.0, 0.0};

    Vector result;
    geometry.ShapeFunctionsValues(result, ksi);

    // Note that the shape function values are evaluated per nodal pair!
    Vector expected_result{3};
    expected_result <<= -0.125, 0.375, 0.75;
    KRATOS_EXPECT_VECTOR_NEAR(result, expected_result, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectAllLocalGradientsAtPosition_ForFourNodedGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateFourNodedLineInterfaceGeometry();
    const auto ksi      = array_1d<double, 3>{0.5, 0.0, 0.0};

    Matrix result;
    geometry.ShapeFunctionsLocalGradients(result, ksi);

    Matrix expected_result(2, 1);
    expected_result <<= -0.5, 0.5;
    KRATOS_EXPECT_MATRIX_NEAR(result, expected_result, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectAllLocalGradientsAtPosition_ForSixNodedGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateSixNodedLineInterfaceGeometry();
    const auto ksi      = array_1d<double, 3>{0.5, 0.0, 0.0};

    Matrix result;
    geometry.ShapeFunctionsLocalGradients(result, ksi);

    Matrix expected_result(3, 1);
    expected_result <<= 0.0, 1.0, -1.0;
    KRATOS_EXPECT_MATRIX_NEAR(result, expected_result, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectJacobian_ForSixNodedGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateSixNodedLineInterfaceGeometry();
    const auto ksi      = array_1d<double, 3>{0.5, 0.0, 0.0};

    Matrix expected_result(2, 1);
    expected_result <<= 3.0, -0.1;

    Matrix result;
    geometry.Jacobian(result, ksi);

    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(result, expected_result, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectJacobian_ForFourNodedGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateFourNodedLineInterfaceGeometry();
    const auto ksi      = array_1d<double, 3>{0.5, 0.0, 0.0};

    Matrix expected_result(2, 1);
    expected_result <<= 3.25, 0.0;

    Matrix result;
    geometry.Jacobian(result, ksi);

    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(result, expected_result, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectDeterminantOfJacobian_ForFourNodedGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateFourNodedLineInterfaceGeometry();
    const auto ksi      = array_1d<double, 3>{0.5, 0.0, 0.0};

    const auto result = geometry.DeterminantOfJacobian(ksi);

    KRATOS_EXPECT_RELATIVE_NEAR(result, 3.25, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectDeterminantOfJacobian_ForSixNodedGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateSixNodedLineInterfaceGeometry();
    const auto ksi      = array_1d<double, 3>{0.5, 0.0, 0.0};

    const auto result = geometry.DeterminantOfJacobian(ksi);

    KRATOS_EXPECT_RELATIVE_NEAR(result, std::sqrt(9.01), 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_Throws_WhenCallingInverseJacobian, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateFourNodedLineInterfaceGeometry();
    const auto ksi      = array_1d<double, 3>{0.5, 0.0, 0.0};

    Matrix result;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        geometry.InverseOfJacobian(result, ksi),
        "Inverse of Jacobian is not implemented for the line interface geometry")
}

KRATOS_TEST_CASE_IN_SUITE(FourNodedLineInterfaceGeometry_LengthReturnsTheLengthOfUnderlyingLineGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateFourNodedLineInterfaceGeometry();

    constexpr auto expected_length = 6.5;
    KRATOS_EXPECT_RELATIVE_NEAR(geometry.Length(), expected_length, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(FourNodedLineInterfaceGeometry_DomainSizeReturnsTheLengthOfUnderlyingLineGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateFourNodedLineInterfaceGeometry();

    constexpr auto expected_length = 6.5;
    KRATOS_EXPECT_RELATIVE_NEAR(geometry.DomainSize(), expected_length, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(SixNodedLineInterfaceGeometry_LengthReturnsTheLengthOfUnderlyingLineGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateSixNodedLineInterfaceGeometry();

    constexpr auto expected_length = 6.504159;
    KRATOS_EXPECT_RELATIVE_NEAR(geometry.Length(), expected_length, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(SixNodedLineInterfaceGeometry_DomainSizeReturnsTheLengthOfUnderlyingLineGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateSixNodedLineInterfaceGeometry();

    constexpr auto expected_length = 6.504159;
    KRATOS_EXPECT_RELATIVE_NEAR(geometry.DomainSize(), expected_length, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(GlobalCoordinatesAreCorrectlyMappedToLocalCoordinateForFourNodedLineInterfaceGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateFourNodedLineInterfaceGeometry();

    auto       result          = array_1d<double, 3>{};
    const auto global_position = array_1d<double, 3>{4.375, 0.1, 0.0};
    result                     = geometry.PointLocalCoordinates(result, global_position);

    const auto expected_result = array_1d<double, 3>{0.5, 0.0, 0.0};
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(result, expected_result, 1e-6)
}

} // namespace Kratos::Testing