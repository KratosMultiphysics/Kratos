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
        return *(std::any_cast<SubRegistryItemPointerType>(mpValue));
    }

    RegistryItem::SubRegistryItemType& RegistryItem::GetSubRegistryItemMap() const
    {
        KRATOS_ERROR_IF(HasValue()) << "Item " << Name() << " has value and cannot be iterated." << std::endl;
        return *(std::any_cast<SubRegistryItemPointerType>(mpValue));
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
        KRATOS_ERROR_IF(iterator == r_map.end()) << "The RegistryItem " << this->Name() << " does not have an item with name " << rItemName << std::endl;
        return *(iterator->second);
    }

    RegistryItem& RegistryItem::GetItem(std::string const& rItemName)
    {
        SubRegistryItemType& r_map = GetSubRegistryItemMap();
        auto iterator = r_map.find(rItemName);
        KRATOS_ERROR_IF(iterator == r_map.end()) << "The RegistryItem " << this->Name() << " does not have an item with name " << rItemName << std::endl;
        return *(iterator->second);
    }

    void RegistryItem::RemoveItem(std::string const& rItemName)
    {
        SubRegistryItemType& r_map = GetSubRegistryItemMap();
        auto iterator = r_map.find(rItemName);
        KRATOS_ERROR_IF(iterator == r_map.end()) << "The RegistryItem " << this->Name() << " does not have an item with name " << rItemName << std::endl;
        r_map.erase(iterator);
    }

} // namespace Kratos.
