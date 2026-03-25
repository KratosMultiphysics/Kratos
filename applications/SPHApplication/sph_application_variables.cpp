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

#include "sph_application_variables.h"

namespace Kratos
{

KRATOS_CREATE_VARIABLE(std::vector<Element::Pointer>, NEIGHBOURS)
KRATOS_CREATE_VARIABLE(double, SMOOTHING_LENGTH)
KRATOS_CREATE_VARIABLE(Matrix, GRADIENT_CORRECTION)
KRATOS_CREATE_VARIABLE(Vector, VW_DKERNEL)
KRATOS_CREATE_VARIABLE(double, VW_KERNEL)
KRATOS_CREATE_VARIABLE(Vector, INTEGRATION_CORRECTION)


KRATOS_CREATE_VARIABLE(double, VOLUME)
KRATOS_CREATE_VARIABLE(Vector, BOUNDARY_NORMAL_AREA)

KRATOS_CREATE_VARIABLE(Vector, SPH_KERNEL_GRADIENT)
KRATOS_CREATE_VARIABLE(double, SPH_KERNEL);


}
