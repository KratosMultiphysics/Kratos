// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Gennady Markelov
//

#include "custom_utilities/check_utilities.h"
#include "includes/checks.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CheckUtilities_CheckDomainSize, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr std::size_t id = 1;
    PointerVector<Node>   nodes;
    nodes.push_back(make_intrusive<Node>(0, 0.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    const auto line_with_coincident_nodes = Line2D2<Node>(nodes);
    nodes.push_back(make_intrusive<Node>(2, 0.0, 0.0, 0.0));
    const auto triangle_with_coincident_nodes = Triangle2D3<Node>(nodes);
    nodes.push_back(make_intrusive<Node>(2, 0.0, 0.0, 0.0));
    const auto tetra_with_coincident_nodes = Tetrahedra3D4<Node>(nodes);

    // Act and Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        CheckUtilities::CheckDomainSize(line_with_coincident_nodes.DomainSize(), id, "Length"),
        "Length (0) is smaller than 1e-15 for element 1")

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        CheckUtilities::CheckDomainSize(triangle_with_coincident_nodes.DomainSize(), id),
        "DomainSize (0) is smaller than 1e-15 for element 1")

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        CheckUtilities::CheckDomainSize(tetra_with_coincident_nodes.DomainSize(), id),
        "DomainSize (0) is smaller than 1e-15 for element 1")
}

} // namespace Kratos::Testing
