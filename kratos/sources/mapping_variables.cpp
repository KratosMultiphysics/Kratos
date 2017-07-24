//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    msandre
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/mapping_variables.h"
#include "includes/kernel.h"

namespace Kratos
{
    KRATOS_CREATE_VARIABLE( boost::shared_ptr<ConditionMap>, CONTACT_MAPS)

void KratosApplication::RegisterMappingVariables()
{
    KRATOS_REGISTER_VARIABLE(CONTACT_MAPS)
}

}  // namespace Kratos.
