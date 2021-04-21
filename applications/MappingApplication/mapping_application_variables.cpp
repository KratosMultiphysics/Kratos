//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

#include "mapping_application_variables.h"

namespace Kratos
{
    KRATOS_CREATE_VARIABLE(int, INTERFACE_EQUATION_ID)
    KRATOS_CREATE_VARIABLE(int, PAIRING_STATUS)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(CURRENT_COORDINATES)
    KRATOS_CREATE_VARIABLE(bool, IS_PROJECTED_LOCAL_SYSTEM)
    KRATOS_CREATE_VARIABLE(bool, IS_DUAL_MORTAR)
}
