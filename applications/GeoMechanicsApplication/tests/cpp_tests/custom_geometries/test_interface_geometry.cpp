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

namespace
{

using namespace Kratos;

LineInterfaceGeometry CreateLineInterfaceGeometry()
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 5.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, -1.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 7.0, 0.0, 0.0));
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

    const auto new_geometry = geometry.Create(nodes);

    KRATOS_EXPECT_NE(new_geometry, nullptr);
    KRATOS_EXPECT_NE(dynamic_cast<const LineInterfaceGeometry*>(new_geometry.get()), nullptr);
    KRATOS_EXPECT_EQ(new_geometry->PointsNumber(), 2);
    KRATOS_EXPECT_EQ(new_geometry->Id(), 0);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_CreateWithId_CreatesNewInstanceOfCorrectTypeAndId,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto          geometry = LineInterfaceGeometry();
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 0.0, 0.0, 0.0));

    const auto new_geometry = geometry.Create(1, nodes);

    KRATOS_EXPECT_NE(new_geometry, nullptr);
    KRATOS_EXPECT_NE(dynamic_cast<const LineInterfaceGeometry*>(new_geometry.get()), nullptr);
    KRATOS_EXPECT_EQ(new_geometry->PointsNumber(), 2);
    KRATOS_EXPECT_EQ(new_geometry->Id(), 1);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectShapeFunctionValuesInNodes_ForTwoLines,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateLineInterfaceGeometry();

    const auto ksi_1 = array_1d<double, 3>{-1.0, 0.0, 0.0};
    const auto ksi_2 = array_1d<double, 3>{1.0, 0.0, 0.0};

    // First pair of nodes
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(0, ksi_1), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(0, ksi_2), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(2, ksi_1), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(2, ksi_2), 0.0);

    // Second pair of nodes
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(1, ksi_1), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(1, ksi_2), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(3, ksi_1), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(3, ksi_2), 1.0);
}

} // namespace Kratos::Testing