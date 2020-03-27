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
    for(auto const& key_variable_pair : KratosComponents<VariableData>::GetComponents()){
        auto variable = key_variable_pair.second;
        if(variable->IsNotComponent()){
            if(KratosComponents<VariableData>::Has(variable->Name() + "_X") && 
               KratosComponents<VariableData>::Has(variable->Name() + "_Y") && 
               KratosComponents<VariableData>::Has(variable->Name() + "_Z")){
                   
                auto const& variable_x = KratosComponents<VariableData>::Get(variable->Name() + "_X");
                auto const& variable_y = KratosComponents<VariableData>::Get(variable->Name() + "_Y");
                auto const& variable_z = KratosComponents<VariableData>::Get(variable->Name() + "_Z");
                if(variable_x.IsComponent() && variable_y.IsComponent() && variable_z.IsComponent()){
                    KRATOS_CHECK_EQUAL(variable_x.Key() + 1, variable_y.Key()) << " for " << variable_x << " and " << variable_y << std::endl;
                    KRATOS_CHECK_EQUAL(variable_y.Key() + 1, variable_z.Key()) << " for " << variable_y << " and " << variable_z << std::endl;
                }          
           }
        }
    }
  }

  KRATOS_TEST_CASE_IN_SUITE(VariableComponent, KratosCoreFastSuite) {
    KRATOS_CHECK_EQUAL(DISPLACEMENT_X.GetSourceVariable(), DISPLACEMENT);
    array_1d<double, 3> displacement = ZeroVector(3);
    displacement[0] = 1.2;
    displacement[1] = 2.3;
    displacement[2] = 3.4;
    KRATOS_CHECK_EQUAL(DISPLACEMENT_X.GetValue(displacement), 1.2);
    KRATOS_CHECK_EQUAL(DISPLACEMENT_Y.GetValue(displacement), 2.3);
    KRATOS_CHECK_EQUAL(DISPLACEMENT_Z.GetValue(displacement), 3.4);

  }


}
} // namespace Kratos.
