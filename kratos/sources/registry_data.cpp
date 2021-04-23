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
#include "includes/registry_data.h"

namespace Kratos
{

    std::string RegistryData::Info() const{
        return mName + " RegistryData";
    }

    void RegistryData::PrintInfo(std::ostream &rOStream) const{
        rOStream << Info();
    }

    void RegistryData::PrintData(std::ostream &rOStream) const{
        
    }

} // namespace Kratos.
