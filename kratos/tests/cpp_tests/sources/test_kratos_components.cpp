//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/kratos_components.h"


namespace Kratos {
namespace Testing {

namespace {

// dummy-conditions required for testing
class DummyCondition1 : public Condition {};
class DummyCondition2 : public Condition {};

}

KRATOS_TEST_CASE_IN_SUITE(KratosComponentsGetNonExistingElement, KratosCoreFastSuite)
{
    KRATOS_EXPECT_TRUE(KratosComponents<Element>::Has("Element2D2N"));
    KRATOS_EXPECT_FALSE(KratosComponents<Element>::Has("NonExisting2D2N"));

    KRATOS_DEBUG_EXCEPT_EXCEPTION_IS_THROWN(KratosComponents<Element>::Get("NonExisting2D2N"), "Error: The component \"NonExisting2D2N\" is not registered!\nMaybe you need to import the application where it is defined?\nThe following components of this type are registered:");
}

KRATOS_TEST_CASE_IN_SUITE(KratosComponentsGetNonExistingVariable, KratosCoreFastSuite)
{
    KRATOS_EXPECT_TRUE(KratosComponents<Variable<double>>::Has("TIME"));
    KRATOS_EXPECT_FALSE(KratosComponents<Variable<double>>::Has("NON_EXISTING_VARIABLE_NAME"));

    KRATOS_DEBUG_EXCEPT_EXCEPTION_IS_THROWN(KratosComponents<Variable<double>>::Get("NON_EXISTING_VARIABLE_NAME"), "Error: The component \"NON_EXISTING_VARIABLE_NAME\" is not registered!\nMaybe you need to import the application where it is defined?\nThe following components of this type are registered:");
}

KRATOS_TEST_CASE_IN_SUITE(KratosComponentsAddDifferentObjectsSameName, KratosCoreFastSuite)
{
    DummyCondition1 dummy_1;
    DummyCondition1 dummy_1_1;

    DummyCondition2 dummy_2;

    // registering the first object - OK
    KratosComponents<Condition>::Add("dummy_1", dummy_1);

    // registering the same object with the same object but with a different name - OK
    KratosComponents<Condition>::Add("dummy_1111", dummy_1);

    // registering different object but of same type - OK
    KratosComponents<Condition>::Add("dummy_1", dummy_1_1);

    // registering the a different object with the name - NOT OK, this is UNDEFINED BEHAVIOR, we don't know what we get when we query the name
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(KratosComponents<Condition>::Add("dummy_1", dummy_2), "Error: An object of different type was already registered with name \"dummy_1\"");

    // Clean up after ourselves
    KratosComponents<Condition>::Remove("dummy_1");
    KratosComponents<Condition>::Remove("dummy_1111");
}

KRATOS_TEST_CASE_IN_SUITE(KratosComponentsRemove, KratosCoreFastSuite)
{
    DummyCondition1 remove_dummy;
    std::string registered_name("RemoveTestDummy");

    // First we add the condition
    KratosComponents<Condition>::Add(registered_name, remove_dummy);
    KRATOS_EXPECT_TRUE(KratosComponents<Condition>::Has(registered_name));

    // Then we remove it
    KratosComponents<Condition>::Remove(registered_name);
    KRATOS_EXPECT_FALSE(KratosComponents<Condition>::Has(registered_name));

    // We can still register another object with the same name
    DummyCondition1 another_dummy;
    KratosComponents<Condition>::Add(registered_name, another_dummy);
    KRATOS_EXPECT_TRUE(KratosComponents<Condition>::Has(registered_name));

    // If we try to remove things that do not exist, we get an error
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        KratosComponents<Condition>::Remove("WrongName"),
        "Error: Trying to remove inexistent component \"WrongName\"."
    );

    // Clean up after ourselves
    KratosComponents<Condition>::Remove(registered_name);
}

}
}  // namespace Kratos.
