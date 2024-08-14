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

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(ConstructInterfaceGeometry, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = LineInterfaceGeometry();
}

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
}

} // namespace Kratos::Testing