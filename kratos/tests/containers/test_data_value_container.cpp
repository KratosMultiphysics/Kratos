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
}
} // namespace Kratos.
