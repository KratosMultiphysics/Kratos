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

#include "custom_elements/U_Pw_small_strain_element.hpp"
#include "custom_utilities/extrapolation_utilities.h"
#include "geometries/triangle_2d_10.h"
#include "includes/checks.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

#include <boost/numeric/ublas/assignment.hpp>

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(ExtrapolationUtilities_CalculateNodalStresses, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    PointerVector<Node> nodes;
    nodes.push_back(make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(2, 1.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(3, 1.0, 1.0, 0.0));
    const auto p_geometry   = std::make_shared<Triangle2D3<Node>>(nodes);
    const auto p_properties = std::make_shared<Properties>();
    auto       element      = UPwSmallStrainElement<2, 3>(1, p_geometry, p_properties, nullptr);

    Vector cauchy_stress_vector(4);
    cauchy_stress_vector <<= 1000.0, 2000.0, 3000.0, 4000.0;
    std::vector<Vector> cauchy_stress_vectors;
    cauchy_stress_vectors.push_back(cauchy_stress_vector);
    cauchy_stress_vector += Vector(4, 1000.0);
    cauchy_stress_vectors.push_back(cauchy_stress_vector);

    std::vector<std::size_t> node_ids = {1, 2, 3};

    // Act and Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(ExtrapolationUtilities::CalculateNodalStresses(
                                          node_ids, element.GetGeometry(), element.GetIntegrationMethod(),
                                          cauchy_stress_vectors, element.Id()),
                                      "An extrapolation matrix size 3 is not equal to given "
                                      "stress vectors size 2 for element Id 1");

    cauchy_stress_vector += Vector(4, 1000.0);
    cauchy_stress_vectors.push_back(cauchy_stress_vector);

    // Act
    auto nodal_stresses = ExtrapolationUtilities::CalculateNodalStresses(
        node_ids, element.GetGeometry(), element.GetIntegrationMethod(), cauchy_stress_vectors,
        element.Id());

    // Assert
    KRATOS_EXPECT_EQ(nodal_stresses.size(), node_ids.size());

    std::vector<Vector> expected_nodal_stresses;
    Vector              expected_stress(4);
    expected_stress = cauchy_stress_vectors[0] - Vector(4, 1000.0);
    expected_nodal_stresses.push_back(expected_stress);
    expected_stress = cauchy_stress_vectors[1];
    expected_nodal_stresses.push_back(expected_stress);
    expected_stress = cauchy_stress_vectors[2] + Vector(4, 1000.0);
    expected_nodal_stresses.push_back(expected_stress);

    for (auto i = std::size_t{0}; i < nodal_stresses.size(); ++i) {
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(nodal_stresses[i].value(), expected_nodal_stresses[i],
                                           Defaults::relative_tolerance);
    }

    // a reduced number of node ids
    node_ids.pop_back();
    nodal_stresses = ExtrapolationUtilities::CalculateNodalStresses(node_ids, element.GetGeometry(),
                                                                    element.GetIntegrationMethod(),
                                                                    cauchy_stress_vectors, element.Id());
    KRATOS_EXPECT_EQ(nodal_stresses.size(), node_ids.size());
    for (auto i = std::size_t{0}; i < nodal_stresses.size(); ++i) {
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(nodal_stresses[i].value(), expected_nodal_stresses[i],
                                           Defaults::relative_tolerance);
    }

    // a wrong node id
    node_ids       = {1, 4, 3};
    nodal_stresses = ExtrapolationUtilities::CalculateNodalStresses(node_ids, element.GetGeometry(),
                                                                    element.GetIntegrationMethod(),
                                                                    cauchy_stress_vectors, element.Id());
    KRATOS_EXPECT_EQ(nodal_stresses.size(), node_ids.size());
    for (auto i = std::size_t{0}; i < nodal_stresses.size(); ++i) {
        if (nodal_stresses[i].has_value()) {
            KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(
                nodal_stresses[i].value(), expected_nodal_stresses[i], Defaults::relative_tolerance);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(ExtrapolationUtilities_CalculateNodalStressesForTriangle2D6, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    PointerVector<Node> nodes;
    nodes.push_back(make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(2, 1.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(3, 1.0, 1.0, 0.0));
    nodes.push_back(make_intrusive<Node>(4, 0.5, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(5, 1.0, 0.5, 0.0));
    nodes.push_back(make_intrusive<Node>(6, 0.5, 0.5, 0.0));
    const auto p_geometry   = std::make_shared<Triangle2D6<Node>>(nodes);
    const auto p_properties = std::make_shared<Properties>();
    auto       element      = UPwSmallStrainElement<2, 6>(1, p_geometry, p_properties, nullptr);

    Vector cauchy_stress_vector(4);
    cauchy_stress_vector <<= 1000.0, 2000.0, 3000.0, 4000.0;
    std::vector<Vector> cauchy_stress_vectors;
    cauchy_stress_vectors.push_back(cauchy_stress_vector);
    cauchy_stress_vector += Vector(4, 1000.0);
    cauchy_stress_vectors.push_back(cauchy_stress_vector);
    cauchy_stress_vector += Vector(4, 1000.0);
    cauchy_stress_vectors.push_back(cauchy_stress_vector);

    // Act
    std::vector<std::size_t> node_ids       = {1, 2, 4};
    auto                     nodal_stresses = ExtrapolationUtilities::CalculateNodalStresses(
        node_ids, element.GetGeometry(), element.GetIntegrationMethod(), cauchy_stress_vectors,
        element.Id());

    // Assert
    KRATOS_EXPECT_EQ(nodal_stresses.size(), node_ids.size());

    std::vector<Vector> expected_nodal_stresses;
    Vector              expected_stress(4);
    expected_stress = cauchy_stress_vectors[0] - Vector(4, 1000.0);
    expected_nodal_stresses.push_back(expected_stress);
    expected_stress = cauchy_stress_vectors[1];
    expected_nodal_stresses.push_back(expected_stress);
    expected_stress = cauchy_stress_vectors[0];
    expected_nodal_stresses.push_back(expected_stress);

    for (auto i = std::size_t{0}; i < nodal_stresses.size(); ++i) {
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(nodal_stresses[i].value(), expected_nodal_stresses[i],
                                           Defaults::relative_tolerance);
    }
}
} // namespace Kratos::Testing
