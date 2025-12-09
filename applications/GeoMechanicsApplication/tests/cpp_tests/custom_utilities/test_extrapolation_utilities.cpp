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

#include "custom_elements/U_Pw_small_strain_element.h"
#include "custom_utilities/extrapolation_utilities.h"
#include "custom_utilities/nodal_extrapolator.h"
#include "geometries/triangle_2d_10.h"
#include "includes/checks.h"
#include "test_setup_utilities/element_setup_utilities.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

#include <boost/numeric/ublas/assignment.hpp>

using namespace Kratos;

namespace Kratos::Testing
{
KRATOS_TEST_CASE_IN_SUITE(ExtrapolationUtilities_CalculateNodalVectors, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_element = ElementSetupUtilities::Create2D3NElement();

    Vector cauchy_stress(4);
    cauchy_stress <<= 1000.0, 2000.0, 3000.0, 4000.0;
    const auto          delta_stress = Vector(4, 1000.0);
    std::vector<Vector> cauchy_stress_vectors;
    cauchy_stress_vectors.push_back(cauchy_stress);
    cauchy_stress_vectors.emplace_back(cauchy_stress + delta_stress);

    std::vector<std::size_t> node_ids = {1, 2, 3};

    // Act and Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        (void)ExtrapolationUtilities::CalculateNodalVectors(node_ids, p_element->GetGeometry(),
                                                            p_element->GetIntegrationMethod(),
                                                            cauchy_stress_vectors, p_element->Id()),
        "An extrapolation matrix size 3 is not equal to given "
        "stress vectors size 2 for element Id 1");

    cauchy_stress_vectors.emplace_back(cauchy_stress + 2 * delta_stress);

    // Act
    auto nodal_stresses = ExtrapolationUtilities::CalculateNodalVectors(
        node_ids, p_element->GetGeometry(), p_element->GetIntegrationMethod(),
        cauchy_stress_vectors, p_element->Id());

    // Assert
    KRATOS_EXPECT_EQ(nodal_stresses.size(), node_ids.size());

    std::vector<Vector> expected_nodal_stresses;
    expected_nodal_stresses.emplace_back(cauchy_stress_vectors[0] - delta_stress);
    expected_nodal_stresses.push_back(cauchy_stress_vectors[1]);
    expected_nodal_stresses.emplace_back(cauchy_stress_vectors[2] + delta_stress);

    for (auto i = std::size_t{0}; i < nodal_stresses.size(); ++i) {
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(nodal_stresses[i].value(), expected_nodal_stresses[i],
                                           Defaults::relative_tolerance);
    }

    // a reduced number of node ids
    node_ids.pop_back();
    nodal_stresses = ExtrapolationUtilities::CalculateNodalVectors(
        node_ids, p_element->GetGeometry(), p_element->GetIntegrationMethod(),
        cauchy_stress_vectors, p_element->Id());
    KRATOS_EXPECT_EQ(nodal_stresses.size(), node_ids.size());
    for (auto i = std::size_t{0}; i < nodal_stresses.size(); ++i) {
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(nodal_stresses[i].value(), expected_nodal_stresses[i],
                                           Defaults::relative_tolerance);
    }

    // a wrong node id
    node_ids       = {1, 4, 3};
    nodal_stresses = ExtrapolationUtilities::CalculateNodalVectors(
        node_ids, p_element->GetGeometry(), p_element->GetIntegrationMethod(),
        cauchy_stress_vectors, p_element->Id());
    KRATOS_EXPECT_EQ(nodal_stresses.size(), node_ids.size());
    KRATOS_EXPECT_FALSE(nodal_stresses[1].has_value())
    for (auto i = std::size_t{0}; i < nodal_stresses.size(); ++i) {
        if (nodal_stresses[i].has_value()) {
            KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(
                nodal_stresses[i].value(), expected_nodal_stresses[i], Defaults::relative_tolerance);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(ExtrapolationUtilities_CalculateNodalVectorsForTriangle2D6, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto p_element = ElementSetupUtilities::Create2D6NElement();

    Vector cauchy_stress(4);
    cauchy_stress <<= 1000.0, 2000.0, 3000.0, 4000.0;
    const auto          delta_stress = Vector(4, 1000.0);
    std::vector<Vector> cauchy_stress_vectors;
    cauchy_stress_vectors.push_back(cauchy_stress);
    cauchy_stress_vectors.emplace_back(cauchy_stress + delta_stress);
    cauchy_stress_vectors.emplace_back(cauchy_stress + 2 * delta_stress);

    // Act
    const std::vector<std::size_t> node_ids       = {1, 4, 2};
    const auto                     nodal_stresses = ExtrapolationUtilities::CalculateNodalVectors(
        node_ids, p_element->GetGeometry(), p_element->GetIntegrationMethod(),
        cauchy_stress_vectors, p_element->Id());

    // Assert
    KRATOS_EXPECT_EQ(nodal_stresses.size(), node_ids.size());

    std::vector<Vector> expected_nodal_stresses;
    expected_nodal_stresses.emplace_back(cauchy_stress_vectors[0] - delta_stress);
    expected_nodal_stresses.push_back(cauchy_stress_vectors[0]);
    expected_nodal_stresses.push_back(cauchy_stress_vectors[1]);

    for (auto i = std::size_t{0}; i < nodal_stresses.size(); ++i) {
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(nodal_stresses[i].value(), expected_nodal_stresses[i],
                                           Defaults::relative_tolerance);
    }
}
} // namespace Kratos::Testing
