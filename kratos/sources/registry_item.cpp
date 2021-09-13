//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

// System includes

// External includes

// Project includes
#include "includes/registry_item.h"

namespace Kratos
{

    std::string RegistryItem::Info() const{
        return mName + " RegistryItem";
    }

    void RegistryItem::PrintInfo(std::ostream &rOStream) const{
        rOStream << Info();
    }

    void RegistryItem::PrintData(std::ostream &rOStream) const{
        
    }

    RegistryItem& RegistryItem::AddItem(Kratos::unique_ptr<RegistryItem>&& pNewItem){
        KRATOS_ERROR_IF(this->HasItem(pNewItem->Name())) << "The RegistryItem " << this->Name() << " already has an item with name " << pNewItem->Name() << std::endl;
        auto insert_result = mSubRegistryItem.insert(std::make_pair(pNewItem->Name(), std::move(pNewItem)));
        KRATOS_ERROR_IF_NOT(insert_result.second) << "Error in inserting " << pNewItem->Name() << " in registry item with name " << this->Name() << std::endl;
        return *insert_result.first->second;
    }


} // namespace Kratos.
