//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
//

#if !defined(KRATOS_FLUID_TRANSPORT_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_FLUID_TRANSPORT_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"
#include "fluid_dynamics_application_variables.h"


namespace Kratos
{
KRATOS_DEFINE_VARIABLE( double, PECLET)
KRATOS_DEFINE_VARIABLE( double, THETA)
KRATOS_DEFINE_VARIABLE( double, PHI_THETA)
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( PHI_GRADIENT)
}

#endif	/* KRATOS_FLUID_TRANSPORT_APPLICATION_VARIABLES_H_INCLUDED */
