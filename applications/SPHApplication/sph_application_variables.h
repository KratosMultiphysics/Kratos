//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "containers/variable.h"
#include "includes/element.h"
#include "containers/global_pointers_vector.h"

namespace Kratos
{

KRATOS_DEFINE_APPLICATION_VARIABLE(SPH_APPLICATION, std::vector<Element::Pointer>, NEIGHBOURS)
KRATOS_DEFINE_APPLICATION_VARIABLE(SPH_APPLICATION, double, SMOOTHING_LENGTH)
KRATOS_DEFINE_APPLICATION_VARIABLE(SPH_APPLICATION, Matrix, GRADIENT_CORRECTION)
KRATOS_DEFINE_APPLICATION_VARIABLE(SPH_APPLICATION, Vector, VW_DKERNEL)
KRATOS_DEFINE_APPLICATION_VARIABLE(SPH_APPLICATION, double, VW_KERNEL)

KRATOS_DEFINE_APPLICATION_VARIABLE(SPH_APPLICATION, double, VOLUME)

KRATOS_DEFINE_APPLICATION_VARIABLE(SPH_APPLICATION, Vector, SPH_KERNEL_GRADIENT)
KRATOS_DEFINE_APPLICATION_VARIABLE(SPH_APPLICATION, double, SPH_KERNEL)


}
