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
#include "custom_utilities/nodal_extrapolator.h"
#include "geometries/triangle_2d_10.h"
#include "includes/checks.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

#include <boost/numeric/ublas/assignment.hpp>

using namespace Kratos;

namespace Kratos::Testing
{

class FakeExtrapolator : public NodalExtrapolator
{
public:
    FakeExtrapolator(std::size_t rows, std::size_t cols) : mRows(rows), mCols(cols) {}

    [[nodiscard]] Matrix CalculateElementExtrapolationMatrix(const GeometryType& rGeometry,
                                                             const GeometryData::IntegrationMethod& rIntegrationMethod) const override
    {
        Matrix mat(mRows, mCols);
        mat.clear();
        return mat;
    }

private:
    std::size_t mRows;
    std::size_t mCols;
};

KRATOS_TEST_CASE_IN_SUITE(ExtrapolationUtilities_CalculateExtrapolationMatrixThrowsExceptions,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    PointerVector<Node> nodes;
    nodes.push_back(make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(2, 1.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(3, 1.0, 1.0, 0.0));
    const auto geometry           = Triangle2D3<Node>(nodes);
    const auto integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;
    const auto fake_row           = FakeExtrapolator(2, 3);
    const auto fake_colums        = FakeExtrapolator(3, 2);

    // Act and Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        auto matrix = ExtrapolationUtilities::CalculateExtrapolationMatrix(geometry, integration_method, 0, &fake_row), "A number of extrapolation matrix rows 2 is not equal to a number of nodes 3 for element id 0");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        auto matrix = ExtrapolationUtilities::CalculateExtrapolationMatrix(geometry, integration_method, 0, &fake_colums), "A number of extrapolation matrix columns 2 is not equal to a number of integration points 3 for element id 0");
}

KRATOS_TEST_CASE_IN_SUITE(ExtrapolationUtilities_CalculateNodalStresses, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    PointerVector<Node> nodes;
    nodes.push_back(make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(2, 1.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(3, 1.0, 1.0, 0.0));
    const auto p_geometry = std::make_shared<Triangle2D3<Node>>(nodes);
    const auto element    = UPwSmallStrainElement<2, 3>(1, p_geometry, nullptr);

    Vector cauchy_stress(4);
    cauchy_stress <<= 1000.0, 2000.0, 3000.0, 4000.0;
    const auto          delta_stress = Vector(4, 1000.0);
    std::vector<Vector> cauchy_stress_vectors;
    cauchy_stress_vectors.push_back(cauchy_stress);
    cauchy_stress_vectors.emplace_back(cauchy_stress + delta_stress);

    std::vector<std::size_t> node_ids = {1, 2, 3};

    // Act and Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(auto stresses = ExtrapolationUtilities::CalculateNodalStresses(
                                          node_ids, element.GetGeometry(), element.GetIntegrationMethod(),
                                          cauchy_stress_vectors, element.Id()),
                                      "An extrapolation matrix size 3 is not equal to given "
                                      "stress vectors size 2 for element Id 1");

    cauchy_stress_vectors.emplace_back(cauchy_stress + 2 * delta_stress);

    // Act
    auto nodal_stresses = ExtrapolationUtilities::CalculateNodalStresses(
        node_ids, element.GetGeometry(), element.GetIntegrationMethod(), cauchy_stress_vectors,
        element.Id());

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
    KRATOS_EXPECT_TRUE(nodal_stresses[1].has_value() == false)
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
    const auto p_geometry = std::make_shared<Triangle2D6<Node>>(nodes);
    const auto element    = UPwSmallStrainElement<2, 6>(1, p_geometry, nullptr);

    Vector cauchy_stress(4);
    cauchy_stress <<= 1000.0, 2000.0, 3000.0, 4000.0;
    const auto          delta_stress = Vector(4, 1000.0);
    std::vector<Vector> cauchy_stress_vectors;
    cauchy_stress_vectors.push_back(cauchy_stress);
    cauchy_stress_vectors.emplace_back(cauchy_stress + delta_stress);
    cauchy_stress_vectors.emplace_back(cauchy_stress + 2 * delta_stress);

    // Act
    std::vector<std::size_t> node_ids       = {1, 4, 2};
    auto                     nodal_stresses = ExtrapolationUtilities::CalculateNodalStresses(
        node_ids, element.GetGeometry(), element.GetIntegrationMethod(), cauchy_stress_vectors,
        element.Id());

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
