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
#include "containers/variables_list_data_value_container.h"
#include "includes/variables.h"
#include "testing/testing.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(VariablesListHas, KratosCoreFastSuite) {
    VariablesList variables_list;
    variables_list.Add(NODAL_AREA);
    variables_list.Add(VELOCITY);


    KRATOS_EXPECT_TRUE(variables_list.Has(NODAL_AREA));
    KRATOS_EXPECT_TRUE(variables_list.Has(VELOCITY));
    KRATOS_EXPECT_FALSE(variables_list.Has(DISPLACEMENT));
}

KRATOS_TEST_CASE_IN_SUITE(VariablesListHasComponent, KratosCoreFastSuite) {
    VariablesList variables_list;
    variables_list.Add(DISPLACEMENT);
 
    KRATOS_EXPECT_FALSE(variables_list.Has(VELOCITY_X));
    KRATOS_EXPECT_FALSE(variables_list.Has(VELOCITY_Y));
    KRATOS_EXPECT_FALSE(variables_list.Has(VELOCITY_Z));

    KRATOS_EXPECT_TRUE(variables_list.Has(DISPLACEMENT_X));
    KRATOS_EXPECT_TRUE(variables_list.Has(DISPLACEMENT_Y));
    KRATOS_EXPECT_TRUE(variables_list.Has(DISPLACEMENT_Z));
}


KRATOS_TEST_CASE_IN_SUITE(VariablesListGetDofInfo, KratosCoreFastSuite) {
    VariablesList variables_list;
    auto dof_index = variables_list.AddDof(&DISPLACEMENT_Y, &REACTION_Y);
    KRATOS_EXPECT_EQ(variables_list.GetDofVariable(dof_index), DISPLACEMENT_Y);
    KRATOS_EXPECT_EQ(variables_list.pGetDofReaction(dof_index), &REACTION_Y);

}

}
} // namespace Kratos.
