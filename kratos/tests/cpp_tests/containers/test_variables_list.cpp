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
#include "includes/stream_serializer.h"
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

namespace {

VariablesList::Pointer CreateSavedDofsList()
{
    auto p_list = Kratos::make_intrusive<VariablesList>();
    p_list->Add(DISPLACEMENT);
    p_list->Add(TEMPERATURE);
    p_list->AddDof(&DISPLACEMENT_X, &REACTION_X);
    p_list->AddDof(&DISPLACEMENT_Y);
    p_list->AddDof(&DISPLACEMENT_Z, &REACTION_Z);
    return p_list;
}

void LoadSavedDofsList(VariablesList::Pointer pList)
{
    StreamSerializer serializer;
    serializer.save("List", CreateSavedDofsList());
    serializer.load("List", pList);
}

} // namespace

// Restart adds the variables before loading and the dofs after it, so the archive's dofs are added to a list with the
// same variables. Variables are compared by key, since each shared library has its own instance of a variable.
KRATOS_TEST_CASE_IN_SUITE(VariablesListLoadAddsDofsToListWithVariables, KratosCoreFastSuite) {
    const Variable<double> temperature_instance(TEMPERATURE);
    auto p_list = Kratos::make_intrusive<VariablesList>();
    p_list->Add(DISPLACEMENT);
    p_list->Add(temperature_instance);
    LoadSavedDofsList(p_list);

    KRATOS_EXPECT_EQ(p_list->size(), 2);
    KRATOS_EXPECT_EQ(p_list->GetDofVariable(0), DISPLACEMENT_X);
    KRATOS_EXPECT_EQ(p_list->pGetDofReaction(0), &REACTION_X);
    KRATOS_EXPECT_EQ(p_list->GetDofVariable(1), DISPLACEMENT_Y);
    KRATOS_EXPECT_EQ(p_list->pGetDofReaction(1), nullptr);
    KRATOS_EXPECT_EQ(p_list->GetDofVariable(2), DISPLACEMENT_Z);
    KRATOS_EXPECT_EQ(p_list->pGetDofReaction(2), &REACTION_Z);

    // a missing reaction is taken from the archive, and a dof added after the save is kept
    auto p_list_with_dofs = Kratos::make_intrusive<VariablesList>();
    p_list_with_dofs->Add(DISPLACEMENT);
    p_list_with_dofs->Add(TEMPERATURE);
    p_list_with_dofs->AddDof(&DISPLACEMENT_X);
    p_list_with_dofs->AddDof(&DISPLACEMENT_Y);
    p_list_with_dofs->AddDof(&DISPLACEMENT_Z, &REACTION_Z);
    p_list_with_dofs->AddDof(&TEMPERATURE);
    LoadSavedDofsList(p_list_with_dofs);

    KRATOS_EXPECT_EQ(p_list_with_dofs->pGetDofReaction(0), &REACTION_X);
    KRATOS_EXPECT_EQ(p_list_with_dofs->pGetDofReaction(1), nullptr);
    KRATOS_EXPECT_EQ(p_list_with_dofs->pGetDofReaction(2), &REACTION_Z);
    KRATOS_EXPECT_EQ(p_list_with_dofs->GetDofVariable(3), TEMPERATURE);
}

// A load into a list with variables validates the archive before changing anything.
KRATOS_TEST_CASE_IN_SUITE(VariablesListLoadRejectsMismatchingArchive, KratosCoreFastSuite) {
    auto p_reordered_list = Kratos::make_intrusive<VariablesList>();
    p_reordered_list->Add(TEMPERATURE);
    p_reordered_list->Add(DISPLACEMENT);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(LoadSavedDofsList(p_reordered_list), "The archive does not match this VariablesList");

    auto p_reordered_dofs_list = Kratos::make_intrusive<VariablesList>();
    p_reordered_dofs_list->Add(DISPLACEMENT);
    p_reordered_dofs_list->Add(TEMPERATURE);
    p_reordered_dofs_list->AddDof(&DISPLACEMENT_Y);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(LoadSavedDofsList(p_reordered_dofs_list), "start with the same dofs, in the same order");

    // the first dof would get a reaction if the archive were applied while checking it
    auto p_other_reaction_list = Kratos::make_intrusive<VariablesList>();
    p_other_reaction_list->Add(DISPLACEMENT);
    p_other_reaction_list->Add(TEMPERATURE);
    p_other_reaction_list->AddDof(&DISPLACEMENT_X);
    p_other_reaction_list->AddDof(&DISPLACEMENT_Y);
    p_other_reaction_list->AddDof(&DISPLACEMENT_Z, &REACTION_Y);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(LoadSavedDofsList(p_other_reaction_list), "The archive gives the dof DISPLACEMENT_Z the reaction REACTION_Z, but this VariablesList has REACTION_Y.");
    KRATOS_EXPECT_EQ(p_other_reaction_list->pGetDofReaction(0), nullptr);
    KRATOS_EXPECT_EQ(p_other_reaction_list->pGetDofReaction(2), &REACTION_Y);
}

}
} // namespace Kratos.
