//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala
//  Collaborators:   Vicente Mataix
//

// System includes

// External includes

// Project includes
#include "includes/master_slave_constraint.h"
#include "includes/kratos_flags.h"

namespace Kratos {

bool MasterSlaveConstraint::IsActive() const 
{
    return IsDefined(ACTIVE) ? Is(ACTIVE) : true;
}

}  // namespace Kratos.