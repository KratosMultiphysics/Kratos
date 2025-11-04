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
#include "custom_utilities/extrapolation_utilities.hpp"
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
    Model               model;
    auto&               r_model_part = model.CreateModelPart("Main");
    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(3, 1.0, 1.0, 0.0));
    const auto p_geometry   = std::make_shared<Triangle2D3<Node>>(nodes);
    const auto p_properties = std::make_shared<Properties>();
    auto       element      = UPwSmallStrainElement<2, 3>(1, p_geometry, p_properties, nullptr);

    Vector cauchy_stress_vector(4);
    cauchy_stress_vector <<= 1000.0, 2000.0, 3000.0, 4000.0;
    std::vector<Vector> cauchy_stress_vectors;
    cauchy_stress_vectors.push_back(cauchy_stress_vector);
    cauchy_stress_vector += Vector(4, 1000.0);
    cauchy_stress_vectors.push_back(cauchy_stress_vector);
    cauchy_stress_vector += Vector(4, 1000.0);
    cauchy_stress_vectors.push_back(cauchy_stress_vector);

    // Act
    const std::vector<std::size_t> node_ids = {1, 2, 3};
    const auto nodal_stresses               = ExtrapolationUtilities::CalculateNodalStresses<3>(
        node_ids, element.GetGeometry(), element.GetIntegrationMethod(), cauchy_stress_vectors);

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
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(nodal_stresses[i], expected_nodal_stresses[i],
                                           Defaults::relative_tolerance);
    }
}

} // namespace Kratos::Testing
