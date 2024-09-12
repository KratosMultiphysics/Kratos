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

#include "custom_geometries/line_interface_geometry.h"
#include "custom_utilities/element_utilities.hpp"
#include "geometries/line_2d_2.h"
#include "includes/node.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

#include <boost/numeric/ublas/assignment.hpp>

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
    const auto geometry = LineInterfaceGeometry<Line2D2<Node>>{1, nodes};

    const Geo::IntegrationPointVectorType integration_points{{-1.0, 0.0, 0.0, 1.0}, {1.0, 0.0, 0.0, 1.0}};
    const auto shape_function_values =
        GeoElementUtilities::EvaluateShapeFunctionsAtIntegrationPoints(integration_points, geometry);

    Vector shape_function_values_1(2);
    shape_function_values_1 <<= 1.0, 0.0;
    Vector shape_function_values_2(2);
    shape_function_values_2 <<= 0.0, 1.0;
    const std::vector<Vector> expected_shape_function_values{shape_function_values_1, shape_function_values_2};

    KRATOS_EXPECT_EQ(expected_shape_function_values.size(), shape_function_values.size());
    for (std::size_t i = 0; i < expected_shape_function_values.size(); ++i) {
        KRATOS_CHECK_VECTOR_NEAR(expected_shape_function_values[i], shape_function_values[i], 1e-6);
    }
}

} // namespace Kratos::Testing
