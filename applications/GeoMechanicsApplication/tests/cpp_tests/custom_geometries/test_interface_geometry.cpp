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
    nodes.push_back(Kratos::make_intrusive<Node>(3, -1.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 7.0, 0.0, 0.0));
    return {1, nodes};
}

LineInterfaceGeometry CreateSixNodedLineInterfaceGeometry()
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 5.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, 2.5, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, -1.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(5, 7.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(6, 3.0, 0.0, 0.0));
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

    // First pair of nodes
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(0, ksi_start), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(0, ksi_end), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(2, ksi_start), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(2, ksi_end), 0.0);

    // Second pair of nodes
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(1, ksi_start), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(1, ksi_end), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(3, ksi_start), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(3, ksi_end), 1.0);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectShapeFunctionValuesInNodes_ForSixNodedGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateSixNodedLineInterfaceGeometry();

    const auto ksi_start  = array_1d<double, 3>{-1.0, 0.0, 0.0};
    const auto ksi_end    = array_1d<double, 3>{1.0, 0.0, 0.0};
    const auto ksi_middle = array_1d<double, 3>{0.0, 0.0, 0.0};

    // First pair of nodes
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(0, ksi_start), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(0, ksi_end), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(0, ksi_middle), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(3, ksi_start), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(3, ksi_end), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(3, ksi_middle), 0.0);

    // Second pair of nodes
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(1, ksi_start), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(1, ksi_end), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(1, ksi_middle), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(4, ksi_start), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(4, ksi_end), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(4, ksi_middle), 0.0);

    // Third pair of nodes
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(2, ksi_start), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(2, ksi_end), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(2, ksi_middle), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(5, ksi_start), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(5, ksi_end), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(5, ksi_middle), 1.0);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectAllShapeFunctionValuesAtPosition_ForFourNodedGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateFourNodedLineInterfaceGeometry();

    const auto ksi = array_1d<double, 3>{0.5, 0.0, 0.0};

    Vector result;
    geometry.ShapeFunctionsValues(result, ksi);

    Vector expected_result{4};
    expected_result <<= 0.25, 0.75, 0.25, 0.75;
    KRATOS_EXPECT_VECTOR_NEAR(result, expected_result, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectAllShapeFunctionValuesAtPosition_ForSixNodedGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateSixNodedLineInterfaceGeometry();

    const auto ksi = array_1d<double, 3>{0.5, 0.0, 0.0};

    Vector result;
    geometry.ShapeFunctionsValues(result, ksi);

    Vector expected_result{6};
    expected_result <<= -0.125, 0.375, 0.75, -0.125, 0.375, 0.75;
    KRATOS_EXPECT_VECTOR_NEAR(result, expected_result, 1e-6)
}

} // namespace Kratos::Testing