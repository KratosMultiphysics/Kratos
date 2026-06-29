//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Anne van de Graaf
//

#include "containers/model.h"
#include "includes/expect.h"
#include "testing/testing.h"
#include "utilities/read_materials_utility.h"

#include <string>

using namespace std::string_literals;

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(ReadMaterialsUtilityRaisesAnErrorWhenAnEmptyStringIsGiven, KratosCoreFastSuiteWithoutKernel) {
    auto model = Model{};
    const auto test_properties = ""s;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN((ReadMaterialsUtility{test_properties, model}), "attempting to parse an empty input; check that your input string or stream contains the expected JSON")
}

KRATOS_TEST_CASE_IN_SUITE(ReadMaterialsUtilityRaisesAnErrorWhenAnEmptyJsonObjectIsGiven, KratosCoreFastSuiteWithoutKernel) {
    auto model = Model{};
    const auto test_properties = "{}"s;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN((ReadMaterialsUtility{test_properties, model}), "Getting a value that does not exist. entry string : properties")
}

KRATOS_TEST_CASE_IN_SUITE(ReadMaterialsUtilityDoesNotRaiseAnErrorWhenTheProperyListIsEmpty, KratosCoreFastSuiteWithoutKernel) {
    auto model = Model{};
    const auto test_properties = R"({
        "properties": []
    })"s;
    EXPECT_NO_THROW((ReadMaterialsUtility{test_properties, model}));
}

KRATOS_TEST_CASE_IN_SUITE(ReadMaterialsUtilityRaisesAnErrorWhenNeitherAModelPartNameNorANameListIsGiven, KratosCoreFastSuiteWithoutKernel) {
    auto model = Model{};
    const auto test_properties = R"({
        "properties": [
            {
                "properties_id": 1
            }
        ]
    })"s;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN((ReadMaterialsUtility{test_properties, model}), "Getting a value that does not exist. entry string : model_part_name")
}

KRATOS_TEST_CASE_IN_SUITE(ReadMaterialsUtilityRaisesAnErrorWhenAnEmptyModelPartNameListIsGiven, KratosCoreFastSuiteWithoutKernel) {
    auto model = Model{};
    const auto test_properties = R"({
        "properties": [
            {
                "properties_id": 1,
                "model_part_name_list": []
            }
        ]
    })"s;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN((ReadMaterialsUtility{test_properties, model}), "At least one model part name must be provided for property 1")
}

KRATOS_TEST_CASE_IN_SUITE(ReadMaterialsUtilityRaisesAnErrorWhenBothModelPartNameAndNameListAreGiven, KratosCoreFastSuiteWithoutKernel) {
    auto model = Model{};
    const auto test_properties = R"({
        "properties": [
            {
                "properties_id": 1,
                "model_part_name": "Foo",
                "model_part_name_list": ["Bar"]
            }
        ]
    })"s;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN((ReadMaterialsUtility{test_properties, model}), "Property 1 provides 'model_part_name' as well as 'model_part_name_list'.")
}

KRATOS_TEST_CASE_IN_SUITE(ReadMaterialsUtilityCreatesAPropertyAndAssignsItToTheModelPart, KratosCoreFastSuiteWithoutKernel) {
    auto model = Model{};
    auto& r_model_part_foo = model.CreateModelPart("Foo"s);
    const auto test_properties = R"({
        "properties": [
            {
                "properties_id": 1,
                "model_part_name": "Foo",
                "Material": {}
            }
        ]
    })"s;

    auto utility = ReadMaterialsUtility{model};
    utility.ReadMaterials(test_properties);

    KRATOS_EXPECT_TRUE(r_model_part_foo.HasProperties(1))
}

KRATOS_TEST_CASE_IN_SUITE(ReadMaterialsUtilityCanCreatePropertiesSharedByModelParts, KratosCoreFastSuiteWithoutKernel) {
    auto model = Model{};
    auto& r_model_part_foo = model.CreateModelPart("Foo"s);
    auto& r_model_part_bar = model.CreateModelPart("Bar"s);
    const auto test_properties = R"({
        "properties": [
            {
                "properties_id": 1,
                "model_part_name_list": ["Foo", "Bar"],
                "Material": {}
            }
        ]
    })"s;

    auto utility = ReadMaterialsUtility{model};
    utility.ReadMaterials(test_properties);

    KRATOS_EXPECT_TRUE(r_model_part_foo.HasProperties(1))
    KRATOS_EXPECT_TRUE(r_model_part_bar.HasProperties(1))
}

KRATOS_TEST_CASE_IN_SUITE(ReadMaterialsUtilityCanCreatePropertiesPerModelPart, KratosCoreFastSuiteWithoutKernel) {
    auto model = Model{};
    auto& r_model_part_foo = model.CreateModelPart("Foo"s);
    auto& r_model_part_bar = model.CreateModelPart("Bar"s);
    const auto test_properties = R"({
        "properties": [
            {
                "properties_id": 1,
                "model_part_name": "Foo",
                "Material": {}
            },
            {
                "properties_id": 2,
                "model_part_name": "Bar",
                "Material": {}
            }
        ]
    })"s;

    auto utility = ReadMaterialsUtility{model};
    utility.ReadMaterials(test_properties);

    KRATOS_EXPECT_TRUE(r_model_part_foo.HasProperties(1))
    KRATOS_EXPECT_TRUE(r_model_part_bar.HasProperties(2))
    KRATOS_EXPECT_FALSE(r_model_part_foo.HasProperties(2))
    KRATOS_EXPECT_FALSE(r_model_part_bar.HasProperties(1))
}

KRATOS_TEST_CASE_IN_SUITE(ReadMaterialsUtilityAssignsPropertiesOfModelPartToElement, KratosCoreFastSuiteWithoutKernel) {
    auto model = Model{};
    auto& r_model_part_foo = model.CreateModelPart("Foo"s);
    constexpr auto id_42 = Element::IndexType{42};
    auto p_element_42 = make_intrusive<Element>(id_42);
    r_model_part_foo.AddElement(p_element_42);
    const auto test_properties = R"({
        "properties": [
            {
                "properties_id": 1,
                "model_part_name": "Foo",
                "Material": {}
            }
        ]
    })"s;

    auto utility = ReadMaterialsUtility{model};
    utility.ReadMaterials(test_properties);

    ASSERT_TRUE(r_model_part_foo.GetElement(id_42).HasProperties());
    KRATOS_EXPECT_EQ(r_model_part_foo.GetElement(id_42).GetProperties().Id(), 1);
}

KRATOS_TEST_CASE_IN_SUITE(ReadMaterialsUtilityAssignsSharedPropertiesOfModelPartsToElements, KratosCoreFastSuiteWithoutKernel) {
    auto model = Model{};

    auto& r_model_part_foo = model.CreateModelPart("Foo"s);
    constexpr auto id_42 = Element::IndexType{42};
    auto p_element_42 = make_intrusive<Element>(id_42);
    r_model_part_foo.AddElement(p_element_42);

    auto& r_model_part_bar = model.CreateModelPart("Bar"s);
    constexpr auto id_43 = Element::IndexType{43};
    auto p_element_43 = make_intrusive<Element>(id_43);
    r_model_part_bar.AddElement(p_element_43);

    const auto test_properties = R"({
        "properties": [
            {
                "properties_id": 1,
                "model_part_name_list": ["Foo", "Bar"],
                "Material": {}
            }
        ]
    })"s;

    auto utility = ReadMaterialsUtility{model};
    utility.ReadMaterials(test_properties);

    ASSERT_TRUE(r_model_part_foo.GetElement(id_42).HasProperties());
    KRATOS_EXPECT_EQ(r_model_part_foo.GetElement(id_42).GetProperties().Id(), 1);
    ASSERT_TRUE(r_model_part_bar.GetElement(id_43).HasProperties());
    KRATOS_EXPECT_EQ(r_model_part_bar.GetElement(id_43).GetProperties().Id(), 1);
}

KRATOS_TEST_CASE_IN_SUITE(ReadMaterialsUtilityAssignsPropertiesOfModelPartToCondition, KratosCoreFastSuiteWithoutKernel) {
    auto model = Model{};
    auto& r_model_part_foo = model.CreateModelPart("Foo"s);
    constexpr auto id_42 = Condition::IndexType{42};
    auto p_condition_42 = make_intrusive<Condition>(id_42);
    r_model_part_foo.AddCondition(p_condition_42);
    const auto test_properties = R"({
        "properties": [
            {
                "properties_id": 1,
                "model_part_name": "Foo",
                "Material": {}
            }
        ]
    })"s;

    auto utility = ReadMaterialsUtility{model};
    utility.ReadMaterials(test_properties);

    ASSERT_TRUE(r_model_part_foo.GetCondition(id_42).HasProperties());
    KRATOS_EXPECT_EQ(r_model_part_foo.GetCondition(id_42).GetProperties().Id(), 1);
}

KRATOS_TEST_CASE_IN_SUITE(ReadMaterialsUtilityAssignsSharedPropertiesOfModelPartsToConditions, KratosCoreFastSuiteWithoutKernel) {
    auto model = Model{};

    auto& r_model_part_foo = model.CreateModelPart("Foo"s);
    constexpr auto id_42 = Condition::IndexType{42};
    auto p_condition_42 = make_intrusive<Condition>(id_42);
    r_model_part_foo.AddCondition(p_condition_42);

    auto& r_model_part_bar = model.CreateModelPart("Bar"s);
    constexpr auto id_43 = Element::IndexType{43};
    auto p_condition_43 = make_intrusive<Condition>(id_43);
    r_model_part_bar.AddCondition(p_condition_43);

    const auto test_properties = R"({
        "properties": [
            {
                "properties_id": 1,
                "model_part_name_list": ["Foo", "Bar"],
                "Material": {}
            }
        ]
    })"s;

    auto utility = ReadMaterialsUtility{model};
    utility.ReadMaterials(test_properties);

    ASSERT_TRUE(r_model_part_foo.GetCondition(id_42).HasProperties());
    KRATOS_EXPECT_EQ(r_model_part_foo.GetCondition(id_42).GetProperties().Id(), 1);
    ASSERT_TRUE(r_model_part_bar.GetCondition(id_43).HasProperties());
    KRATOS_EXPECT_EQ(r_model_part_bar.GetCondition(id_43).GetProperties().Id(), 1);
}

}