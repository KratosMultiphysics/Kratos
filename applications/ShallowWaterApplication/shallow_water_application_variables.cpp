//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#include "shallow_water_application_variables.h"

namespace Kratos
{
// Water depth and bathymetry data
KRATOS_CREATE_VARIABLE(int,PATCH_INDEX)    // TODO: Remove variable ??
KRATOS_CREATE_VARIABLE(double,ELEVATION)   // Free surface elevation (eta)
KRATOS_CREATE_VARIABLE(double,BATHYMETRY)  // Bathymetry (H)
KRATOS_CREATE_VARIABLE(double,DEPTH)       // Depth (h=H+eta)
}
