//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "containers/variable.h"

namespace Kratos
{
    KRATOS_DEFINE_APPLICATION_VARIABLE( MAPPING_APPLICATION, int, INTERFACE_EQUATION_ID ) // Has to be int bcs of MPI
    KRATOS_DEFINE_APPLICATION_VARIABLE( MAPPING_APPLICATION, int, PAIRING_STATUS )
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( MAPPING_APPLICATION, CURRENT_COORDINATES )
    KRATOS_DEFINE_APPLICATION_VARIABLE( MAPPING_APPLICATION, bool, IS_PROJECTED_LOCAL_SYSTEM)
    KRATOS_DEFINE_APPLICATION_VARIABLE(MAPPING_APPLICATION, bool, IS_DUAL_MORTAR)
}
