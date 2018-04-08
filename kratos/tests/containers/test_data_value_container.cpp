//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
// kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Vicente Mataix Ferrandiz
//
//

// System includes
#include <unordered_set>

// External includes

// Project includes
#include "containers/data_value_container.h"
#include "includes/variables.h"
#include "testing/testing.h"

namespace Kratos {
namespace Testing {

/**
 * This test checks the proper functioning of the DataValueContainer
 */
KRATOS_TEST_CASE_IN_SUITE(DataValueContainer, KratosCoreFastSuite) {
  DataValueContainer container;
  Vector original_distances(4);
  original_distances[0] = 0.00;
  original_distances[1] = 0.10;
  original_distances[2] = 0.20;
  original_distances[3] = 0.30;
  container.SetValue(ELEMENTAL_DISTANCES, original_distances);
  auto& distances = container.GetValue(ELEMENTAL_DISTANCES);

  for (std::size_t i = 0; i < distances.size(); i++)
    KRATOS_CHECK_EQUAL(distances[i], original_distances[i]);
}

/**
 * This test checks recursively the proper functioning of the DataValueContainer
 */
KRATOS_TEST_CASE_IN_SUITE(DataValueContainerRecursive, KratosCoreFastSuite) {
  DataValueContainer container;
  Vector original_distances(4);
  original_distances[0] = 0.00;
  original_distances[1] = 0.10;
  original_distances[2] = 0.20;
  original_distances[3] = 0.30;

  const std::size_t number_recursions = 1000;
  for (std::size_t i = 0; i < number_recursions; ++i) {
    container.SetValue(NODAL_H, 0.5);
    container.SetValue(ELEMENTAL_DISTANCES, original_distances);
    const auto& nodal_h = container.GetValue(NODAL_H);
    const auto& distances = container.GetValue(ELEMENTAL_DISTANCES);

    KRATOS_CHECK_EQUAL(nodal_h, 0.5);
    for (std::size_t i = 0; i < distances.size(); i++) {
        KRATOS_CHECK_EQUAL(distances[i], original_distances[i]);
    }
  }
}

/**
 * This test checks recursively the proper functioning of the DataValueContainer for components varioables
 */
KRATOS_TEST_CASE_IN_SUITE(DataValueContainerRecursiveComponents, KratosCoreFastSuite) {
  DataValueContainer container;
  array_1d<double, 3> aux_zero_array_value(3, 0.0);

  const std::size_t number_recursions = 10;
  for (std::size_t i = 0; i < number_recursions; ++i) {
    container.SetValue(DISPLACEMENT, aux_zero_array_value);
    container.SetValue(DISPLACEMENT_Y, 1.0);
    container.SetValue(DISPLACEMENT_Z, 2.0);
    const auto& displacement = container.GetValue(DISPLACEMENT);

    KRATOS_CHECK_EQUAL(displacement[0], 0.0);
    KRATOS_CHECK_EQUAL(displacement[1], 1.0);
    KRATOS_CHECK_EQUAL(displacement[2], 2.0);
  }
}
}
} // namespace Kratos.
