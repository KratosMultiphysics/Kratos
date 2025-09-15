//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// System includes

// External includes

// Project includes
#include "includes/geometrical_object.h"
#include "includes/kratos_flags.h"

namespace Kratos {

bool GeometricalObject::IsActive() const 
{
    return IsDefined(ACTIVE) ? Is(ACTIVE) : true;
}

}  // namespace Kratos.
