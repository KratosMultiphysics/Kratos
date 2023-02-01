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

    RegistryItem::SubRegistryItemType::iterator RegistryItem::begin()
    {
        KRATOS_ERROR_IF(HasValue()) << "Item " << Name() << " has value and cannot be iterated." << std::endl;
        return mSubRegistryItem.begin();
    }

    RegistryItem::SubRegistryItemType::const_iterator RegistryItem::cbegin() const
    {
        KRATOS_ERROR_IF(HasValue()) << "Item " << Name() << " has value and cannot be iterated." << std::endl;
        return mSubRegistryItem.cbegin();
    }

    RegistryItem::SubRegistryItemType::iterator RegistryItem::end()
    {
        KRATOS_ERROR_IF(HasValue()) << "Item " << Name() << " has value and cannot be iterated." << std::endl;
        return mSubRegistryItem.end();
    }

    RegistryItem::SubRegistryItemType::const_iterator RegistryItem::cend() const
    {
        KRATOS_ERROR_IF(HasValue()) << "Item " << Name() << " has value and cannot be iterated." << std::endl;
        return mSubRegistryItem.cend();
    }

    RegistryItem::KeyReturnConstIterator RegistryItem::KeyConstBegin() const
    {
        KRATOS_ERROR_IF(HasValue()) << "Item " << Name() << " has value and cannot be iterated." << std::endl;
        return KeyReturnConstIterator(mSubRegistryItem.cbegin());
    }

    RegistryItem::KeyReturnConstIterator RegistryItem::KeyConstEnd() const
    {
        KRATOS_ERROR_IF(HasValue()) << "Item " << Name() << " has value and cannot be iterated." << std::endl;
        return KeyReturnConstIterator(mSubRegistryItem.cend());
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
        for(auto& item : mSubRegistryItem){
            rOStream << *(item.second) << std::endl;
        }
    }

    std::string RegistryItem::ToJson(std::string const& rIndentation) const
    {
        KRATOS_ERROR_IF(HasValue()) << "For storing a value you should use the RegistryValueItem" << std::endl;

        std::stringstream buffer;
        buffer  << rIndentation << "\"" << mName << "\" : {" << std::endl;
        for(auto& item : mSubRegistryItem){
            buffer << item.second->ToJson(rIndentation + "    ");
        }
        buffer << "}" << std::endl;

        return buffer.str();
    }

    RegistryItem const& RegistryItem::GetItem(std::string const& rItemName) const
    {
        auto iterator = mSubRegistryItem.find(rItemName);
        KRATOS_ERROR_IF(iterator == mSubRegistryItem.end()) << "The RegistryItem " << this->Name() << " does not have an item with name " << rItemName << std::endl;
        return *(iterator->second);
    }

    RegistryItem& RegistryItem::GetItem(std::string const& rItemName)
    {
        auto iterator = mSubRegistryItem.find(rItemName);
        KRATOS_ERROR_IF(iterator == mSubRegistryItem.end()) << "The RegistryItem " << this->Name() << " does not have an item with name " << rItemName << std::endl;
        return *(iterator->second);
    }

    void RegistryItem::RemoveItem(std::string const& rItemName)
    {
        auto iterator = mSubRegistryItem.find(rItemName);
        KRATOS_ERROR_IF(iterator == mSubRegistryItem.end()) << "The RegistryItem " << this->Name() << " does not have an item with name " << rItemName << std::endl;
        mSubRegistryItem.erase(iterator);
    }


} // namespace Kratos.
