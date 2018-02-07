//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
//kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//

// System includes
#include <unordered_set>

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/variables.h"

namespace Kratos {
namespace Testing {

// This test is to check the hash function used in variable key, if it fails 
// means that we should change the hash function. Pooyan.
KRATOS_TEST_CASE_IN_SUITE(VariablesKeyUniqueness, KratosCoreFastSuite) {
  std::unordered_set<VariableData::KeyType> registered_keys;
  for (auto &i : KratosComponents<VariableData>::GetComponents()) {
    auto hash = i.second->Key();
    if (registered_keys.find(hash) == registered_keys.end()) {
      registered_keys.insert(hash);
    } else {
      KRATOS_ERROR << "Duplicated variables key founded for variable" << *(i.second) << std::endl;
    }
  }
}

KRATOS_TEST_CASE_IN_SUITE(VariablesKeyOrder, KratosCoreFastSuite) {
    KRATOS_CHECK_EQUAL(VELOCITY_X.Key() + 1, VELOCITY_Y.Key()) << VELOCITY_X << " , " << VELOCITY_Y;
    KRATOS_CHECK_EQUAL(VELOCITY_Y.Key() + 1, VELOCITY_Z.Key()); 

    KRATOS_CHECK_EQUAL(DISPLACEMENT_X.Key() + 1, DISPLACEMENT_Y.Key());
    KRATOS_CHECK_EQUAL(DISPLACEMENT_Y.Key() + 1, DISPLACEMENT_Z.Key()); 
}


}
} // namespace Kratos.
