//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:   Philipp Bucher
//
//

// System includes

// External includes

// Project includes
#include "includes/properties.h"
#include "utilities/properties_accessor.h"

namespace Kratos
{

Properties::AccessorPointerType& Properties::GetAccessor()
{
    if (!mpAccessor)
        mpAccessor = Kratos::make_shared<PropertiesAccessor>(this);
    return mpAccessor;
}

}  // namespace Kratos.
