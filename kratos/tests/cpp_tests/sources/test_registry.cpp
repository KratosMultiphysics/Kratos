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
#include "includes/checks.h"
#include "includes/registry.h"

namespace Kratos {

namespace Testing {

namespace
{
    void CheckItemIteration(
        const RegistryItem& rRegistryItem,
        const std::vector<std::string>& rRegistryItemKeys)
    {
        std::size_t n_found = 0;
        for (auto it = rRegistryItem.cbegin(); it != rRegistryItem.cend(); ++it) {
            const auto it_found = std::find(rRegistryItemKeys.begin(), rRegistryItemKeys.end(), it->first);
            const bool is_found = it_found != rRegistryItemKeys.end() ? true : false;
            if (is_found) {n_found++;}
            KRATOS_CHECK(is_found);
        }
        KRATOS_CHECK_EQUAL(n_found, rRegistryItemKeys.size());
    }

    void CheckKeyIteration(
        RegistryItem& rRegistryItem,
        const std::vector<std::string>& rRegistryItemKeys)
    {
        std::size_t n_found = 0;
        for (auto it = rRegistryItem.KeyConstBegin(); it != rRegistryItem.KeyConstEnd(); ++it) {
            const auto it_found = std::find(rRegistryItemKeys.begin(), rRegistryItemKeys.end(), *it);
            const bool is_found = it_found != rRegistryItemKeys.end() ? true : false;
            if (is_found) {n_found++;}
        }
        KRATOS_CHECK_EQUAL(n_found, rRegistryItemKeys.size());
    }
}

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

    auto& sub_item = registry_item.GetItem("sub_item");
    KRATOS_CHECK_STRING_EQUAL(sub_item.Name(),"sub_item");
    KRATOS_CHECK_STRING_EQUAL(registry_item.ToJson(), "{\n\"items\": {\n\"sub_item\": {}\n}\n}");

    registry_item.RemoveItem("sub_item");
    KRATOS_CHECK_IS_FALSE(registry_item.HasItems());
    KRATOS_CHECK_IS_FALSE(registry_item.HasItem("sub_item"));
}

KRATOS_TEST_CASE_IN_SUITE(RegistryValue, KratosCoreFastSuite)
{
    RegistryItem empty_value_item("empty_value_item");
    KRATOS_CHECK_STRING_EQUAL(empty_value_item.Name(),"empty_value_item");
    KRATOS_CHECK_IS_FALSE(empty_value_item.HasValue());
    KRATOS_CHECK_IS_FALSE(empty_value_item.HasItems());
    KRATOS_CHECK_IS_FALSE(empty_value_item.HasItem("test"));

    double value = 3.14;
    RegistryItem value_registry_item("value_item", value);
    KRATOS_CHECK_STRING_EQUAL(value_registry_item.Name(),"value_item");
    KRATOS_CHECK(value_registry_item.HasValue());
    KRATOS_CHECK_IS_FALSE(value_registry_item.HasItems());
    KRATOS_CHECK_IS_FALSE(value_registry_item.HasItem("test"));

    KRATOS_CHECK_DOUBLE_EQUAL(value_registry_item.GetValue<double>(), 3.14);

    std::string value_item_json = R"({
"value_item": "3.14"
})";
    KRATOS_CHECK_STRING_EQUAL(value_registry_item.ToJson(), value_item_json);

    RegistryItem registry_item("items");
    registry_item.AddItem<double>("sub_value_item", value);
    KRATOS_CHECK_IS_FALSE(registry_item.HasValue());
    KRATOS_CHECK(registry_item.HasItems());
    KRATOS_CHECK_IS_FALSE(registry_item.HasItem("test"));
    KRATOS_CHECK(registry_item.HasItem("sub_value_item"));

    KRATOS_CHECK_STRING_EQUAL(registry_item.ToJson(), "{\n\"items\": {\n\"sub_value_item\": \"3.14\"\n}\n}");

    auto& sub_item = registry_item.GetItem("sub_value_item");
    KRATOS_CHECK_STRING_EQUAL(sub_item.Name(),"sub_value_item");
    KRATOS_CHECK(sub_item.HasValue());
    registry_item.RemoveItem("sub_value_item");
    KRATOS_CHECK_IS_FALSE(registry_item.HasItems());
    KRATOS_CHECK_IS_FALSE(registry_item.HasItem("sub_value_item"));
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

KRATOS_TEST_CASE_IN_SUITE(RegistryIteration, KratosCoreFastSuite)
{
    Registry::AddItem<RegistryItem>("item_in_root_1");
    Registry::AddItem<RegistryItem>("item_in_root_1.subitem_1");
    Registry::AddItem<RegistryItem>("item_in_root_1.subitem_2");
    Registry::AddItem<RegistryItem>("item_in_root_1.subitem_1.subsubitem_1");
    Registry::AddItem<RegistryItem>("item_in_root_1.subitem_1.subsubitem_2");
    Registry::AddItem<RegistryItem>("item_in_root_1.subitem_1.subsubitem_3");
    Registry::AddItem<RegistryItem>("item_in_root_1.subitem_1.subsubitem_2.subsubsubitem_1");
    Registry::AddItem<RegistryItem>("item_in_root_1.subitem_1.subsubitem_2.subsubsubitem_2");
    Registry::AddItem<RegistryItem>("item_in_root_1.subitem_1.subsubitem_2.subsubsubitem_3");
    Registry::AddItem<RegistryItem>("item_in_root_1.subitem_1.subsubitem_2.subsubsubitem_4");

    // Check subitems iteration
    std::vector<std::string> sub_keys = {"subitem_1","subitem_2"};
    CheckItemIteration(Registry::GetItem("item_in_root_1"), sub_keys);

    // Check subsubitems iteration
    std::vector<std::string> sub_sub_keys = {"subsubitem_1","subsubitem_2","subsubitem_3"};
    CheckItemIteration(Registry::GetItem("item_in_root_1.subitem_1"), sub_sub_keys);

    // Check subsubsubitems iteration
    std::vector<std::string> sub_sub_sub_keys = {"subsubsubitem_1","subsubsubitem_2","subsubsubitem_3","subsubsubitem_4"};
    CheckItemIteration(Registry::GetItem("item_in_root_1.subitem_1.subsubitem_2"), sub_sub_sub_keys);

    // Remove the auxiliary testing items
    Registry::RemoveItem("item_in_root_1");
}

KRATOS_TEST_CASE_IN_SUITE(RegistryKeyIteration, KratosCoreFastSuite)
{
    Registry::AddItem<RegistryItem>("item_in_root_1");
    Registry::AddItem<RegistryItem>("item_in_root_1.subitem_1");
    Registry::AddItem<RegistryItem>("item_in_root_1.subitem_2");
    Registry::AddItem<RegistryItem>("item_in_root_1.subitem_1.subsubitem_1");
    Registry::AddItem<RegistryItem>("item_in_root_1.subitem_1.subsubitem_2");
    Registry::AddItem<RegistryItem>("item_in_root_1.subitem_1.subsubitem_3");
    Registry::AddItem<RegistryItem>("item_in_root_1.subitem_1.subsubitem_2.subsubsubitem_1");
    Registry::AddItem<RegistryItem>("item_in_root_1.subitem_1.subsubitem_2.subsubsubitem_2");
    Registry::AddItem<RegistryItem>("item_in_root_1.subitem_1.subsubitem_2.subsubsubitem_3");
    Registry::AddItem<RegistryItem>("item_in_root_1.subitem_1.subsubitem_2.subsubsubitem_4");

    // Check subitems iteration
    std::vector<std::string> sub_keys = {"subitem_1","subitem_2"};
    CheckKeyIteration(Registry::GetItem("item_in_root_1"), sub_keys);

    // Check subsubitems iteration
    std::vector<std::string> sub_sub_keys = {"subsubitem_1","subsubitem_2","subsubitem_3"};
    CheckKeyIteration(Registry::GetItem("item_in_root_1.subitem_1"), sub_sub_keys);

    // Check subsubsubitems iteration
    std::vector<std::string> sub_sub_sub_keys = {"subsubsubitem_1","subsubsubitem_2","subsubsubitem_3","subsubsubitem_4"};
    CheckKeyIteration(Registry::GetItem("item_in_root_1.subitem_1.subsubitem_2"), sub_sub_sub_keys);

    // Remove the auxiliary testing items
    Registry::RemoveItem("item_in_root_1");
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

KRATOS_TEST_CASE_IN_SUITE(RegistrySomeRegisteredVariables, KratosCoreFastSuite){
     KRATOS_CHECK(Registry::HasItem("variables.all.TEMPERATURE"));
     KRATOS_CHECK(Registry::HasItem("variables.all.VELOCITY"));
     KRATOS_CHECK(Registry::HasItem("variables.all.DISPLACEMENT_Z"));
}


}
}  // namespace Kratos.
