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
#include "includes/cfd_variables.h"
#include "fluid_dynamics_application_variables.h"


namespace Kratos
{
KRATOS_DEFINE_APPLICATION_VARIABLE(FLUID_TRANSPORT_APPLICATION, double, PECLET)
KRATOS_DEFINE_APPLICATION_VARIABLE(FLUID_TRANSPORT_APPLICATION, double, THETA)
KRATOS_DEFINE_APPLICATION_VARIABLE(FLUID_TRANSPORT_APPLICATION, double, PHI_THETA)
KRATOS_DEFINE_APPLICATION_VARIABLE(FLUID_TRANSPORT_APPLICATION, double, NODAL_ANALYTIC_SOLUTION)
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(FLUID_TRANSPORT_APPLICATION, PHI_GRADIENT)
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(FLUID_TRANSPORT_APPLICATION, NODAL_PHI_GRADIENT )

}

#endif	/* KRATOS_FLUID_TRANSPORT_APPLICATION_VARIABLES_H_INCLUDED */
