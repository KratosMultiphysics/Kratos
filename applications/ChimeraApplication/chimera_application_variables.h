//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Authors:        Aditya Ghantasala, https://github.com/adityaghantasala
// 					Navaneeth K Narayanan
//					Rishith Ellath Meethal
//

#if !defined(KRATOS_CHIMERA_APPLICATION_VARIABLES_H_INCLUDED)
#define KRATOS_CHIMERA_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "containers/variable.h"

namespace Kratos
{

KRATOS_DEFINE_APPLICATION_VARIABLE(CHIMERA_APPLICATION,double, CHIMERA_DISTANCE);
KRATOS_DEFINE_APPLICATION_VARIABLE(CHIMERA_APPLICATION,double, ROTATIONAL_ANGLE);
KRATOS_DEFINE_APPLICATION_VARIABLE(CHIMERA_APPLICATION,double, ROTATIONAL_VELOCITY);
KRATOS_DEFINE_APPLICATION_VARIABLE(CHIMERA_APPLICATION,bool, CHIMERA_INTERNAL_BOUNDARY);
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(CHIMERA_APPLICATION, ROTATION_MESH_DISPLACEMENT);
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(CHIMERA_APPLICATION, ROTATION_MESH_VELOCITY);
}

#endif /* KRATOS_CHIMERA_APPLICATION_VARIABLES_H_INCLUDED */
