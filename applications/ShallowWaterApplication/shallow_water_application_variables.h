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


#if !defined(KRATOS_SHALLOW_WATER_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_SHALLOW_WATER_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"

namespace Kratos
{
// Water depth and bathymetry data
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, int,PATCH_INDEX) // TODO: remove variable?
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double,ETA)         // Water surface elevation (eta)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double,BATHYMETRY)  // Bathymetry (H)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double,TOTAL_DEPTH) // Total depth (h=H+eta)
}

#endif	/* KRATOS_FLUID_DYNAMICS_APPLICATION_VARIABLES */
