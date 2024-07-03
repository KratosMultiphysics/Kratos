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
//                   Ruben Zorrilla
//

// System includes
#include <mutex>

// External includes

// Project includes
#include "includes/registry.h"

namespace Kratos
{

namespace
{
    std::once_flag flag_once;
}

    RegistryItem* Registry::mspRootRegistryItem = nullptr;

    RegistryItem& Registry::GetItem(std::string const& rItemFullName)
    {
        const std::lock_guard<LockObject> scope_lock(ParallelUtilities::GetGlobalLock());

        auto item_path = StringUtilities::SplitStringByDelimiter(rItemFullName, '.');
        KRATOS_ERROR_IF(item_path.empty()) << "The item full name is empty" << std::endl;


        RegistryItem* p_current_item = &GetRootRegistryItem();

        for(std::size_t i = 0 ; i < item_path.size() ; i++){
            auto& r_item_name = item_path[i];
            if(p_current_item->HasItem(r_item_name)){
                p_current_item = &p_current_item->GetItem(r_item_name);
            }
            else{
                KRATOS_ERROR << "The item \"" << rItemFullName << "\" is not found in the registry. The item \"" << p_current_item->Name() << "\" does not have \"" << r_item_name << "\"" << std::endl;
            }
        }

        return *p_current_item;
    }

    void Registry::RemoveItem(std::string const& rItemFullName)
    {
        const std::lock_guard<LockObject> scope_lock(ParallelUtilities::GetGlobalLock());

        auto item_path = StringUtilities::SplitStringByDelimiter(rItemFullName, '.');
        KRATOS_ERROR_IF(item_path.empty()) << "The item full name is empty" << std::endl;

        RegistryItem* p_current_item = &GetRootRegistryItem();

        for(std::size_t i = 0 ; i < item_path.size() - 1 ; i++){
            auto& r_item_name = item_path[i];
            if(p_current_item->HasItem(r_item_name)){
                p_current_item = &p_current_item->GetItem(r_item_name);
            }
            else{
                KRATOS_ERROR << "The item \"" << rItemFullName << "\" is not found in the registry. The item \"" << p_current_item->Name() << "\" does not have \"" << r_item_name << "\"" << std::endl;
            }
        }

        auto& r_item_name = item_path.back();
        if(p_current_item->HasItem(r_item_name)){
            p_current_item->RemoveItem(r_item_name);
        }
        else{
            KRATOS_ERROR << "The item \"" << rItemFullName << "\" is not found in the registry. The item \"" << p_current_item->Name() << "\" does not have \"" << r_item_name << "\"" << std::endl;
        }
    }

    std::size_t Registry::size()
    {
        return mspRootRegistryItem->size();
    }

    bool Registry::HasItem(std::string const& rItemFullName)
    {
        const std::lock_guard<LockObject> scope_lock(ParallelUtilities::GetGlobalLock());

        auto item_path = StringUtilities::SplitStringByDelimiter(rItemFullName, '.');
        KRATOS_ERROR_IF(item_path.empty()) << "The item full name is empty" << std::endl;

        RegistryItem* p_current_item = &GetRootRegistryItem();

        for(std::size_t i = 0 ; i < item_path.size() ; i++){
            auto& r_item_name = item_path[i];
            if(p_current_item->HasItem(r_item_name)){
                p_current_item = &p_current_item->GetItem(r_item_name);
            }
            else{
               return false;
            }
        }
        return true;
    }

    bool Registry::HasValue(std::string const& rItemFullName)
    {
        return GetItem(rItemFullName).HasValue();
    }

    bool Registry::HasItems(std::string const& rItemFullName)
    {
        return GetItem(rItemFullName).HasItems();
    }

    std::string Registry::Info() const
    {
        return "Registry";
    }

    void Registry::PrintInfo(std::ostream &rOStream) const
    {
        rOStream << Info();
    }

    void Registry::PrintData(std::ostream &rOStream) const
    {}

    std::string Registry::ToJson(std::string const& Indentation) const
    {
        return GetRootRegistryItem().ToJson(Indentation);
    }

    RegistryItem& Registry::GetRootRegistryItem()
    {
        if (!mspRootRegistryItem) {
            std::call_once(flag_once, [](){
                static RegistryItem root_item("Registry");
                mspRootRegistryItem = &root_item;
            });
        }

        return *mspRootRegistryItem;
    }

} // namespace Kratos.
