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




} // namespace Kratos.
