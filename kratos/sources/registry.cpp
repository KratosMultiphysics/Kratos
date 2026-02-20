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
#include <regex>
#include <memory>

// External includes

// Project includes
#include "includes/registry.h"

// Demangle OS includes
#if defined(__GNUG__) || defined(__clang__)
    #include <cxxabi.h>
#endif

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
            } else {
                NotFoundError(rItemFullName, r_item_name, p_current_item);
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
            } else {
                NotFoundError(rItemFullName, r_item_name, p_current_item);
            }
        }

        auto& r_item_name = item_path.back();
        if(p_current_item->HasItem(r_item_name)){
            p_current_item->RemoveItem(r_item_name);
        } else {
            NotFoundError(rItemFullName, r_item_name, p_current_item);
        }
    }

    void Registry::SetCurrentSource(std::string const & rCurrentSource)
    {
        // If context key not present, create it
        if (Registry::HasItem("CurrentContext")){
            Registry::RemoveItem("CurrentContext");
        }

        // It is needed to create a std::string explicitly copying the '"CurrentContext"+rCurrentSource' to avoid casting problems
        // involing std::any_cast to a reference type which key references a string that may not be alive when invoked.
        std::string context_key = std::string("CurrentContext." + rCurrentSource);

        Registry::AddItem<RegistryItem>(context_key);
    }

    std::string Registry::GetCurrentSource()
    {
        // If context key not present, create it
        if (!Registry::HasItem("CurrentContext")){
            Registry::AddItem<RegistryItem>("CurrentContext.KratosMultiphysics");
        }

        return Registry::GetItem("CurrentContext").begin()->first;
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

    std::string Registry::GetRegistryName(std::string d_name) {
        // Remove namespace 
        static const std::regex ns_re(R"(\b[\w_]+::)");
        d_name = std::regex_replace(d_name, ns_re, "");

        // Remove literal suffix
        static const std::regex suffix_re(R"(\b(\d+)(ull|ul|ll|l|u|f|L)\b)");
        d_name = std::regex_replace(d_name, suffix_re, "$1");

        // Strip all spaces
        static const std::regex space_re(R"(\s+)");
        d_name = std::regex_replace(d_name, space_re, "");

        // Strip the trailing * from the usage of pointer types (Some demanglers leave that, some other do not)
        d_name = std::regex_replace(d_name, std::regex("(>?)\\*$"), "$1");

        return d_name;
    }

#if defined(__GNUG__) || defined(__clang__)
    std::string Registry::GetDemangledName(const char* name) {
        int status = 0;
        std::unique_ptr<char, void(*)(void*)> res {
            abi::__cxa_demangle(name, nullptr, nullptr, &status),
            std::free
        };
        return (status == 0) ? Registry::GetRegistryName(res.get()) : name;
    }
#else
    std::string Registry::GetDemangledName(const char* name) { return Registry::GetRegistryName(name); }
#endif


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

    void Registry::NotFoundError(
        const std::string& rFullName,
        const std::string& rItemName,
        RegistryItem* pCurrentItem
        )
    {
        const std::vector<std::string> available_list = pCurrentItem->GetSubItemAvailableList();
        std::stringstream error_message_buffer;
        error_message_buffer << "The item \"" << rFullName << "\" is not found in the registry. The item \"" << pCurrentItem->Name() << "\" does not have \"" << rItemName << "\". The available objects are: \n";
        for (std::string const& item : available_list) {
            error_message_buffer << "\t\t" << item << "\n";
        }
        KRATOS_ERROR << error_message_buffer.str() << std::endl;
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
