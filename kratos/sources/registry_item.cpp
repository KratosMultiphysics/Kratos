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

// External includes

// Project includes
#include "includes/registry_item.h"

namespace Kratos
{
    RegistryItem::SubRegistryItemType& RegistryItem::GetSubRegistryItemMap()
    {
        KRATOS_ERROR_IF(HasValue()) << "Item " << Name() << " has value and cannot be iterated." << std::endl;
#if defined (__GNUC__) && __GNUC__ <= 8 && __GNUC_MINOR__ <= 3
        return *(boost::any_cast<SubRegistryItemPointerType>(mpValue));
#else
        return *(std::any_cast<SubRegistryItemPointerType>(mpValue));
#endif
    }

    RegistryItem::SubRegistryItemType& RegistryItem::GetSubRegistryItemMap() const
    {
        KRATOS_ERROR_IF(HasValue()) << "Item " << Name() << " has value and cannot be iterated." << std::endl;
#if defined (__GNUC__) && __GNUC__ <= 8 && __GNUC_MINOR__ <= 3
        return *(boost::any_cast<SubRegistryItemPointerType>(mpValue));
#else
        return *(std::any_cast<SubRegistryItemPointerType>(mpValue));
#endif
    }

    RegistryItem::SubRegistryItemType::iterator RegistryItem::begin()
    {
        return GetSubRegistryItemMap().begin();
    }

    RegistryItem::SubRegistryItemType::const_iterator RegistryItem::cbegin() const
    {
        return GetSubRegistryItemMap().cbegin();
    }

    RegistryItem::SubRegistryItemType::iterator RegistryItem::end()
    {
        return GetSubRegistryItemMap().end();
    }

    RegistryItem::SubRegistryItemType::const_iterator RegistryItem::cend() const
    {
        return GetSubRegistryItemMap().cend();
    }

    RegistryItem::KeyReturnConstIterator RegistryItem::KeyConstBegin() const
    {
        return KeyReturnConstIterator(GetSubRegistryItemMap().cbegin());
    }

    RegistryItem::KeyReturnConstIterator RegistryItem::KeyConstEnd() const
    {
        return KeyReturnConstIterator(GetSubRegistryItemMap().cend());
    }

    std::string RegistryItem::Info() const
    {
        return mName + " RegistryItem ";
    }

    void RegistryItem::PrintInfo(std::ostream &rOStream) const
    {
        rOStream << Info();
    }

    void RegistryItem::PrintData(std::ostream &rOStream) const
    {
        if (HasValue()) {
            rOStream << this->GetValueString();
        } else {
            for(auto& item : GetSubRegistryItemMap()){
                rOStream << item.second->GetValueString() << std::endl;
            }
        }
    }

    std::vector<std::string> RegistryItem::GetSubItemAvailableList() const
    {
        std::vector<std::string> available_list;
        auto it_item_begin = this->cbegin();
        auto it_item_end = this->cend();
        available_list.reserve(std::distance(it_item_begin, it_item_end));
        for (auto it_item = it_item_begin; it_item != it_item_end; ++it_item) {
            available_list.push_back((it_item->second)->Name());
        }
        return available_list;
    }

    void RegistryItem::NotFoundError(const std::string& rItemName) const
    {
        const std::vector<std::string> available_list = GetSubItemAvailableList();
        std::stringstream available_list_str;
        available_list_str << "The RegistryItem " << this->Name() << " does not have an item with name " << rItemName << ". The available objects are: \n";
        for (std::string const& item : available_list) {
            available_list_str << "\t\t" << item << "\n";
        }
        KRATOS_ERROR << available_list_str.str() << std::endl;
    }

    std::string RegistryItem::GetValueString() const
    {
        return (this->*(this->mGetValueStringMethod))();
    }

    std::string RegistryItem::ToJson(std::string const& rTabSpacing, const std::size_t Level) const
    {
        std::string tabbing;
        for (std::size_t i = 0; i < Level + 1; ++i) {
            tabbing += rTabSpacing;
        }

        std::stringstream buffer;

        if (Level == 0)
            buffer << "{" << std::endl;

        if (HasValue()) {
            buffer << tabbing << "\"" << mName << "\": \""
                   << this->GetValueString() << "\"";
        } else {
            buffer << tabbing << "\"" << mName << "\": {";

            for (auto& r_item : GetSubRegistryItemMap()) {
                buffer << std::endl;
                buffer << r_item.second->ToJson(rTabSpacing, Level + 1);
                buffer << ",";
            }

            if (HasItems()) {
                buffer.seekp(-1, std::ios_base::end);
                buffer << std::endl << tabbing;
            }

            buffer << "}";
        }

        if (Level == 0)
            buffer << std::endl << "}";

        return buffer.str();
    }

    bool RegistryItem::HasValue() const
    {
        return (mpValue.type() != typeid(SubRegistryItemPointerType));
    }

    bool RegistryItem::HasItem(std::string const& rItemName) const
    {
        if (!HasValue()) {
            SubRegistryItemType& r_map = GetSubRegistryItemMap();
            return (r_map.find(rItemName) != r_map.end());
        } else {
            return false;
        }
    }

    bool RegistryItem::HasItems() const
    {
        return (!HasValue() && !GetSubRegistryItemMap().empty());
    }

    std::size_t RegistryItem::size()
    {
        return GetSubRegistryItemMap().size();
    }

    RegistryItem const& RegistryItem::GetItem(std::string const& rItemName) const
    {
        SubRegistryItemType& r_map = GetSubRegistryItemMap();
        auto iterator = r_map.find(rItemName);
        if (iterator == r_map.end()) {
            NotFoundError(rItemName);
        }
        return *(iterator->second);
    }

    RegistryItem& RegistryItem::GetItem(std::string const& rItemName)
    {
        SubRegistryItemType& r_map = GetSubRegistryItemMap();
        auto iterator = r_map.find(rItemName);
        if (iterator == r_map.end()) {
            NotFoundError(rItemName);
        }
        return *(iterator->second);
    }

    void RegistryItem::RemoveItem(std::string const& rItemName)
    {
        SubRegistryItemType& r_map = GetSubRegistryItemMap();
        auto iterator = r_map.find(rItemName);
        if (iterator == r_map.end()) {
            NotFoundError(rItemName);
        }
        r_map.erase(iterator);
    }

} // namespace Kratos.
