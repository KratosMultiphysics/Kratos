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
#include "includes/kratos_components.h"
#include "includes/condition.h"


namespace Kratos {
namespace Testing {

namespace {

// dummy-conditions required for testing
class DummyCondition1 : public Condition {};
class DummyCondition2 : public Condition {};

}

KRATOS_TEST_CASE_IN_SUITE(KratosComponentsGetNonExistingElement, KratosCoreFastSuite)
{
    KRATOS_CHECK(KratosComponents<Element>::Has("Element2D2N"));
    KRATOS_CHECK_IS_FALSE(KratosComponents<Element>::Has("NonExisting2D2N"));

    KRATOS_DEBUG_CHECK_EXCEPTION_IS_THROWN(KratosComponents<Element>::Get("NonExisting2D2N"), "Error: The component \"NonExisting2D2N\" is not registered!\nMaybe you need to import the application where it is defined?\nThe following components of this type are registered:");
}

KRATOS_TEST_CASE_IN_SUITE(KratosComponentsGetNonExistingVariable, KratosCoreFastSuite)
{
    KRATOS_CHECK(KratosComponents<Variable<double>>::Has("TIME"));
    KRATOS_CHECK_IS_FALSE(KratosComponents<Variable<double>>::Has("NON_EXISTING_VARIABLE_NAME"));

    KRATOS_DEBUG_CHECK_EXCEPTION_IS_THROWN(KratosComponents<Variable<double>>::Get("NON_EXISTING_VARIABLE_NAME"), "Error: The component \"NON_EXISTING_VARIABLE_NAME\" is not registered!\nMaybe you need to import the application where it is defined?\nThe following components of this type are registered:");
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
    KRATOS_CHECK_EXCEPTION_IS_THROWN(KratosComponents<Condition>::Add("dummy_1", dummy_2), "Error: An object of different type was already registered with name \"dummy_1\"");
}

}
}  // namespace Kratos.
