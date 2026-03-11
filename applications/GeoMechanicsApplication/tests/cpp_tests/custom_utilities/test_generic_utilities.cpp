// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Mohamed Nabi
//                   Wijtze Pieter Kikstra
//                   Richard Faasse
//                   Gennady Markelov
//

#include "containers/pointer_vector.h"
#include "custom_utilities/generic_utilities.hpp"
#include "custom_utilities/ublas_utilities.h"
#include "geometries/geometry.h"
#include "includes/expect.h"
#include "includes/node.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.h"

using namespace Kratos;

namespace Kratos::Testing
{

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CheckVectorPermutation)
{
    // Arrange
    const auto vector  = UblasUtilities::CreateVector({1.0, 2.0, 3.0, 4.0, 5.0});
    const auto indices = std::vector<int>{4, 3, 2, 1, 0};

    // Act & assert
    const auto expected_result = UblasUtilities::CreateVector({5.0, 4.0, 3.0, 2.0, 1.0});
    KRATOS_EXPECT_VECTOR_EQ(GenericUtilities::PermutedVector(vector, indices), expected_result);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CheckMatrixPermutation)
{
    // Arrange
    const auto matrix = UblasUtilities::CreateMatrix(
        {{1.0, 0.0, 0.0, 0.0}, {0.0, 1.0, 0.0, 0.0}, {0.0, 0.0, 1.0, 0.0}, {0.0, 0.0, 0.0, 1.0}});
    const auto indices = std::vector<int>{3, 2, 1, 0};

    // Act & assert
    const auto expected_result = UblasUtilities::CreateMatrix(
        {{0.0, 0.0, 0.0, 1.0}, {0.0, 0.0, 1.0, 0.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 0.0, 0.0, 0.0}});
    KRATOS_EXPECT_MATRIX_EQ(GenericUtilities::MatrixWithPermutedColumns(matrix, indices), expected_result);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, GetIdsFromEntityContents_ReturnsEmptyListForEmptyGeometry)
{
    const auto node_ids = GenericUtilities::GetIdsFromEntityContents(Geometry<Node>{});

    KRATOS_EXPECT_TRUE(node_ids.empty())
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, GetIdsFromEntityContents_ReturnsCorrectNodeIds)
{
    // Arange
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(42, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(314, 0.0, 0.0, 0.0));
    const Geometry geometry(1, nodes);

    // Act
    const auto node_ids_vector = GenericUtilities::GetIdsFromEntityContents(geometry);

    std::set<Node::IndexType> node_ids_set;
    GenericUtilities::GetIdsFromEntityContents(geometry, std::inserter(node_ids_set, node_ids_set.end()));

    // Assert
    const std::vector<Node::IndexType> expected_ids_vec = {1, 3, 42, 314};
    const std::set<Node::IndexType> expected_ids_set(expected_ids_vec.begin(), expected_ids_vec.end());

    EXPECT_EQ(node_ids_vector, expected_ids_vec);
    EXPECT_EQ(node_ids_set, expected_ids_set);
}

} // namespace Kratos::Testing