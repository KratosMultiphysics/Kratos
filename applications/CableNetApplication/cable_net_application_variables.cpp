//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Klaus B. Sautter
//

#include "cable_net_application_variables.h"

namespace Kratos
{
    KRATOS_CREATE_VARIABLE(double, RING_TENSILE_STIFFNESS)
    KRATOS_CREATE_VARIABLE(double, RING_BENDING_STIFFNESS)
    KRATOS_CREATE_VARIABLE(double, RING_THICKNESS_WIRE)
    KRATOS_CREATE_VARIABLE(double, RING_NR_WIRES)
    KRATOS_CREATE_VARIABLE(double, RING_REFERENCE_CIRCUMFERENCE)
    KRATOS_CREATE_VARIABLE(double, BRAKING_MASS)
    KRATOS_CREATE_VARIABLE(Vector, SPRING_DEFORMATION_EMPIRICAL_POLYNOMIAL)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(WALL_POSITION)
}
