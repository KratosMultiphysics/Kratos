//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//

// System includes
#include <numeric>

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/registry.h"

namespace Kratos {

namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(RegistryItem, KratosCoreFastSuite)
{
    RegistryItem empty_registry_item("empty_item");
    KRATOS_CHECK_STRING_EQUAL(empty_registry_item.Name(),"empty_item");
    KRATOS_CHECK_IS_FALSE(empty_registry_item.HasValue());
    KRATOS_CHECK_IS_FALSE(empty_registry_item.HasItems());
    KRATOS_CHECK_IS_FALSE(empty_registry_item.HasItem("test"));

    double value = 3.14;
    RegistryItem value_registry_item("value_item", value);
    KRATOS_CHECK_STRING_EQUAL(value_registry_item.Name(),"value_item");
    KRATOS_CHECK(value_registry_item.HasValue());
    KRATOS_CHECK_IS_FALSE(value_registry_item.HasItems());
    KRATOS_CHECK_IS_FALSE(value_registry_item.HasItem("test"));

    RegistryItem registry_item("items");
    registry_item.AddItem<RegistryItem>("sub_item");
    KRATOS_CHECK_IS_FALSE(registry_item.HasValue());
    KRATOS_CHECK(registry_item.HasItems());
    KRATOS_CHECK_IS_FALSE(registry_item.HasItem("test"));
    KRATOS_CHECK(registry_item.HasItem("sub_item"));

    std::cout << registry_item << std::endl;

    auto& sub_item = registry_item.GetItem("sub_item");
    KRATOS_CHECK_STRING_EQUAL(sub_item.Name(),"sub_item");
    std::string item_json = R"("items" : {
    "sub_item" : {
}
}
)";
    KRATOS_CHECK_STRING_EQUAL(registry_item.ToJson(), item_json);

    registry_item.RemoveItem("sub_item");
    KRATOS_CHECK_IS_FALSE(registry_item.HasItems());
    KRATOS_CHECK_IS_FALSE(registry_item.HasItem("sub_item"));
}

KRATOS_TEST_CASE_IN_SUITE(RegistryValue, KratosCoreFastSuite)
{
    RegistryValueItem<double> empty_value_item("empty_value_item");
    KRATOS_CHECK_STRING_EQUAL(empty_value_item.Name(),"empty_value_item");
    KRATOS_CHECK_IS_FALSE(empty_value_item.HasValue());
    KRATOS_CHECK_IS_FALSE(empty_value_item.HasItems());
    KRATOS_CHECK_IS_FALSE(empty_value_item.HasItem("test"));

    double value = 3.14;
    RegistryValueItem<double> value_registry_item("value_item", value);
    KRATOS_CHECK_STRING_EQUAL(value_registry_item.Name(),"value_item");
    KRATOS_CHECK(value_registry_item.HasValue());
    KRATOS_CHECK_IS_FALSE(value_registry_item.HasItems());
    KRATOS_CHECK_IS_FALSE(value_registry_item.HasItem("test"));

    std::string value_item_json = R"("value_item" : "3.14"
)";
    KRATOS_CHECK_STRING_EQUAL(value_registry_item.ToJson(), value_item_json);

    RegistryItem registry_item("items");
    registry_item.AddItem<RegistryValueItem<double>>("sub_value_item", value);
    KRATOS_CHECK_IS_FALSE(registry_item.HasValue());
    KRATOS_CHECK(registry_item.HasItems());
    KRATOS_CHECK_IS_FALSE(registry_item.HasItem("test"));
    KRATOS_CHECK(registry_item.HasItem("sub_value_item"));

    std::string item_json = R"("items" : {
    "sub_value_item" : "3.14"
}
)";
    KRATOS_CHECK_STRING_EQUAL(registry_item.ToJson(), item_json);

    auto& sub_item = registry_item.GetItem("sub_value_item");
    KRATOS_CHECK_STRING_EQUAL(sub_item.Name(),"sub_value_item");
    KRATOS_CHECK(sub_item.HasValue());
    registry_item.RemoveItem("sub_value_item");
    KRATOS_CHECK_IS_FALSE(registry_item.HasItems());
    KRATOS_CHECK_IS_FALSE(registry_item.HasItem("sub_value_item"));

    KRATOS_WATCH(registry_item.ToJson());
}


KRATOS_TEST_CASE_IN_SUITE(RegistryAddAndRemove, KratosCoreFastSuite)
{
    KRATOS_CHECK_IS_FALSE(Registry::HasItem("item_in_root"));
    KRATOS_CHECK_IS_FALSE(Registry::HasItem("path.to.the.registry.new_item"));


    Registry::AddItem<RegistryItem>("item_in_root");
    KRATOS_CHECK(Registry::HasItem("item_in_root"));
    KRATOS_CHECK_IS_FALSE(Registry::HasItem("path.to.the.registry.new_item"));
    auto& item_in_root = Registry::GetItem("item_in_root");
    KRATOS_CHECK_STRING_EQUAL(item_in_root.Name(),"item_in_root");


    Registry::AddItem<RegistryItem>("path.to.the.registry.new_item");
    KRATOS_CHECK(Registry::HasItem("item_in_root"));
    KRATOS_CHECK(Registry::HasItem("path.to.the.registry.new_item"));
    auto& new_item = Registry::GetItem("path.to.the.registry.new_item");
    KRATOS_CHECK_STRING_EQUAL(new_item.Name(),"new_item");

    Registry::RemoveItem("item_in_root");
    KRATOS_CHECK_IS_FALSE(Registry::HasItem("item_in_root"));
    KRATOS_CHECK(Registry::HasItem("path.to.the.registry.new_item"));

    Registry::RemoveItem("path.to.the.registry.new_item");
    KRATOS_CHECK_IS_FALSE(Registry::HasItem("item_in_root"));
    KRATOS_CHECK_IS_FALSE(Registry::HasItem("path.to.the.registry.new_item"));
}

KRATOS_TEST_CASE_IN_SUITE(RegistryParallelAddAndRemove, KratosCoreFastSuite)
{
    std::size_t size = 1000;

    IndexPartition<int>(size).for_each(
        [&](int i){
                std::string item_name = "item_" + std::to_string(i);
                KRATOS_CHECK_IS_FALSE(Registry::HasItem(item_name));
            }
        );

    IndexPartition<int>(size).for_each(
        [&](int i){
                std::string item_name = "item_" + std::to_string(i);
                Registry::AddItem<RegistryItem>(item_name);
                KRATOS_CHECK(Registry::HasItem(item_name));
                auto& item = Registry::GetItem(item_name);
                KRATOS_CHECK_STRING_EQUAL(item.Name(),item_name);
            }
        );

    IndexPartition<int>(size).for_each(
        [&](int i){
                std::string item_name = "item_" + std::to_string(i);
                KRATOS_CHECK(Registry::HasItem(item_name));
                auto& item = Registry::GetItem(item_name);
                KRATOS_CHECK_STRING_EQUAL(item.Name(),item_name);
            }
        );

    IndexPartition<int>(size).for_each(
        [&](int i){
                std::string item_name = "item_" + std::to_string(i);
                std::string item_path = std::string("path.to.the.registry.new_item.") + item_name;
                Registry::AddItem<RegistryItem>(item_path);
                KRATOS_CHECK(Registry::HasItem(item_path));
                auto& item = Registry::GetItem(item_path);
                KRATOS_CHECK_STRING_EQUAL(item.Name(),item_name);
            }
        );

    IndexPartition<int>(size).for_each(
        [&](int i){
                std::string item_name = "item_" + std::to_string(i);
                std::string item_path = std::string("path.to.the.registry.new_item.") + item_name;
                KRATOS_CHECK(Registry::HasItem(item_path));
                auto& item = Registry::GetItem(item_path);
                KRATOS_CHECK_STRING_EQUAL(item.Name(),item_name);
            }
        );

    IndexPartition<int>(size).for_each(
        [&](int i){
                std::string item_name = "item_" + std::to_string(i);
                Registry::RemoveItem(item_name);
                KRATOS_CHECK_IS_FALSE(Registry::HasItem(item_name));
            }
        );

    IndexPartition<int>(size).for_each(
        [&](int i){
                std::string item_name = std::string("path.to.the.registry.new_item.item_") + std::to_string(i);
                Registry::RemoveItem(item_name);
                KRATOS_CHECK_IS_FALSE(Registry::HasItem(item_name));
            }
        );

}

}
}  // namespace Kratos.
