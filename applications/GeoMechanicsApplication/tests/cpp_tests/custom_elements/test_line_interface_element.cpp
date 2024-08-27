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

#include "custom_elements/line_interface_element.h"
#include "custom_geometries/line_interface_geometry.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

namespace Kratos::Testing
{

using namespace Kratos;

KRATOS_TEST_CASE_IN_SUITE(LineInterfaceElementIsAnElement, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    LineInterfaceElement element;
    auto                 casted_element = dynamic_cast<Element*>(&element);
    KRATOS_CHECK_NOT_EQUAL(casted_element, nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(LineInterfaceElementCanCreateInstance, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const LineInterfaceElement element;

    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 1.0, 0.0, 0.0));
    const auto geometry = std::make_shared<LineInterfaceGeometry>(nodes);

    auto properties = std::make_shared<Properties>();

    // Act
    auto created_element = element.Create(1, geometry, properties);

    // Assert
    EXPECT_NE(created_element, nullptr);
    EXPECT_EQ(created_element->Id(), 1);
    EXPECT_NE(created_element->pGetGeometry(), nullptr);
    EXPECT_NE(created_element->pGetProperties(), nullptr);
}

} // namespace Kratos::Testing
