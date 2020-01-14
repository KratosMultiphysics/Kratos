//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
// ==============================================================================
//  ChimeraApplication
//
//  License:         BSD License
//                   license: ChimeraApplication/license.txt
//
//  Authors:        Aditya Ghantasala, https://github.com/adityaghantasala
// 					Navaneeth K Narayanan
//					Rishith Ellath Meethal
// ==============================================================================
//

#if !defined(KRATOS_CHIMERA_APPLICATION_VARIABLES_H_INCLUDED)
#define KRATOS_CHIMERA_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"

namespace Kratos
{

KRATOS_DEFINE_VARIABLE(double, CHIMERA_DISTANCE);
KRATOS_DEFINE_VARIABLE(double, ROTATIONAL_ANGLE);
KRATOS_DEFINE_VARIABLE(double, ROTATIONAL_VELOCITY);
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(CHIMERA_APPLICATION, ROTATION_MESH_DISPLACEMENT);
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(CHIMERA_APPLICATION, ROTATION_MESH_VELOCITY);

KRATOS_DEFINE_FLAG( FS_CHIMERA_VELOCITY_CONSTRAINT);
KRATOS_DEFINE_FLAG( FS_CHIMERA_PRESSURE_CONSTRAINT);
KRATOS_DEFINE_FLAG( CHIMERA_INTERNAL_BOUNDARY);
}

#endif /* KRATOS_CHIMERA_APPLICATION_VARIABLES_H_INCLUDED */
