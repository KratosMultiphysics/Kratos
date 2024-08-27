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

namespace
{

using namespace Kratos;

PointerVector<Node> CreateNodes()
{
    PointerVector<Node> result;
    result.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    result.push_back(Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0));
    result.push_back(Kratos::make_intrusive<Node>(3, 0.0, 0.0, 0.0));
    result.push_back(Kratos::make_intrusive<Node>(4, 1.0, 0.0, 0.0));

    return result;
}

} // namespace

namespace Kratos::Testing
{

using namespace Kratos;

KRATOS_TEST_CASE_IN_SUITE(LineInterfaceElementIsAnElement, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    LineInterfaceElement element;
    auto                 casted_element = dynamic_cast<Element*>(&element);
    KRATOS_CHECK_NOT_EQUAL(casted_element, nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(LineInterfaceElementCanCreateInstanceWithGeometryInput, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const LineInterfaceElement element;
    const auto                 geometry   = std::make_shared<LineInterfaceGeometry>(CreateNodes());
    auto                       properties = std::make_shared<Properties>();

    // Act
    auto created_element = element.Create(1, geometry, properties);

    // Assert
    EXPECT_NE(created_element, nullptr);
    EXPECT_EQ(created_element->Id(), 1);
    EXPECT_NE(created_element->pGetGeometry(), nullptr);
    EXPECT_NE(created_element->pGetProperties(), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(LineInterfaceElementCanCreateInstanceWithNodeInput, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto nodes      = CreateNodes();
    auto properties = std::make_shared<Properties>();

    // The source element needs to have a geometry, otherwise the version of the
    // Create method with a node input will fail.
    const auto                 geometry = std::make_shared<LineInterfaceGeometry>(nodes);
    const LineInterfaceElement element(1, geometry, properties);

    // Act
    auto created_element = element.Create(1, nodes, properties);

    // Assert
    EXPECT_NE(created_element, nullptr);
    EXPECT_EQ(created_element->Id(), 1);
    EXPECT_NE(created_element->pGetGeometry(), nullptr);
    EXPECT_NE(created_element->pGetProperties(), nullptr);
}

} // namespace Kratos::Testing
