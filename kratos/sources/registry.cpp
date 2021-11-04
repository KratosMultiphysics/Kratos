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
#include <mutex>

// External includes

// Project includes
#include "includes/registry.h"

namespace Kratos
{

namespace {
    std::once_flag flag_once;
}

    RegistryItem* Registry::mspRootRegistryItem = nullptr;

    std::string Registry::Info() const{
        return "Registry";
    }

    void Registry::PrintInfo(std::ostream &rOStream) const{
        rOStream << Info();
    }

    void Registry::PrintData(std::ostream &rOStream) const{
    }

    std::string Registry::ToJson(std::string const& Indentation) const {
        return const_cast<Registry*>(this)->GetRootRegistryItem().ToJson(Indentation);
    }

    RegistryItem& Registry::GetRootRegistryItem(){
        if (!mspRootRegistryItem) {
            std::call_once(flag_once, [](){
                static RegistryItem root_item("Registry");
                mspRootRegistryItem = &root_item;
            });
        }

    return *mspRootRegistryItem;

    }


} // namespace Kratos.
