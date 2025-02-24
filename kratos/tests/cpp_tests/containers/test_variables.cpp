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
//                   Vicente Mataix Ferrandiz
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

KRATOS_TEST_CASE_IN_SUITE(TestTimeDerivatives, KratosCoreFastSuite)
{
    KRATOS_EXPECT_EQ((DISPLACEMENT.GetTimeDerivative()).Name(), "VELOCITY");
    KRATOS_EXPECT_EQ((DISPLACEMENT_X.GetTimeDerivative()).Name(), "VELOCITY_X");
    KRATOS_EXPECT_EQ(((DISPLACEMENT.GetTimeDerivative()).GetTimeDerivative()).Name(), "ACCELERATION");
    KRATOS_EXPECT_EQ((VELOCITY.GetTimeDerivative()).Name(), "ACCELERATION");
    KRATOS_EXPECT_EQ((ROTATION.GetTimeDerivative()).Name(), "ANGULAR_VELOCITY");
    KRATOS_EXPECT_EQ((ROTATION_X.GetTimeDerivative()).Name(), "ANGULAR_VELOCITY_X");
    KRATOS_EXPECT_EQ(((ROTATION.GetTimeDerivative()).GetTimeDerivative()).Name(), "ANGULAR_ACCELERATION");
    KRATOS_EXPECT_EQ((ANGULAR_VELOCITY.GetTimeDerivative()).Name(), "ANGULAR_ACCELERATION");
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
                    KRATOS_EXPECT_EQ(variable_x.Key() + 1, variable_y.Key()) << " for " << variable_x << " and " << variable_y << std::endl;
                    KRATOS_EXPECT_EQ(variable_y.Key() + 1, variable_z.Key()) << " for " << variable_y << " and " << variable_z << std::endl;
                }
           }
        }
    }
  }

KRATOS_TEST_CASE_IN_SUITE(VariablesRegister, KratosCoreFastSuite) {
    Variable<double> new_var("NEW_DUMMY_VARIABLE");
    new_var.Register();
    KRATOS_EXPECT_TRUE(Registry::HasItem("variables.all.NEW_DUMMY_VARIABLE"));
    KRATOS_EXPECT_TRUE(Registry::HasItem("variables.KratosMultiphysics.NEW_DUMMY_VARIABLE"));

    // Attempting to register a variable with the same name and type does nothing
    Variable<double> duplicate_variable("NEW_DUMMY_VARIABLE");
    duplicate_variable.Register();
    KRATOS_EXPECT_TRUE(Registry::HasItem("variables.all.NEW_DUMMY_VARIABLE"));
    KRATOS_EXPECT_TRUE(Registry::HasItem("variables.KratosMultiphysics.NEW_DUMMY_VARIABLE"));

    // Attempting to register a variable with the same name and different type is an error
    Variable<int> wrong_type_variable("NEW_DUMMY_VARIABLE");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        wrong_type_variable.Register(),
        "Attempting to register NEW_DUMMY_VARIABLE but a variable with the same name and different type already exists"
    )
}


}
} // namespace Kratos.
