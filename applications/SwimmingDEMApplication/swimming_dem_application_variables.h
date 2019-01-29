//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//
//

#ifndef KRATOS_SWIMMING_DEM_APPLICATION_VARIABLES_H
#define	KRATOS_SWIMMING_DEM_APPLICATION_VARIABLES_H

#include "includes/define.h"
#include "includes/variables.h"

namespace Kratos {

    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(SWIMMING_DEM_APPLICATION, AVERAGED_FLUID_VELOCITY)
    KRATOS_DEFINE_APPLICATION_VARIABLE(SWIMMING_DEM_APPLICATION, std::string, SDEM_HYDRODYNAMIC_INTERACTION_LAW_NAME)

}

#endif // KRATOS_SWIMMING_DEM_APPLICATION_VARIABLES_H  defined