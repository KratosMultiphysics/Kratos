//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#if !defined(KRATOS_STABILIZED_CFD_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_STABILIZED_CFD_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/kratos_application.h"

#include "custom_utilities/turbulence_statistics_container.h"

namespace Kratos
{
KRATOS_DEFINE_VARIABLE( double, FIC_BETA )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( DIRECTIONAL_BETA )
KRATOS_DEFINE_VARIABLE( int, RECORDED_STEPS )
KRATOS_DEFINE_VARIABLE( double, MEAN_KINETIC_ENERGY )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( MEAN_VELOCITY )
KRATOS_DEFINE_VARIABLE( double, MEAN_PRESSURE )
KRATOS_DEFINE_VARIABLE( Matrix, VELOCITY_COVARIANCES )
KRATOS_DEFINE_VARIABLE( TurbulenceStatisticsContainer::Pointer, TURBULENCE_STATISTICS )
KRATOS_DEFINE_VARIABLE( double, TRACE_XI )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( DIV_XI )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( MOMENTUM_PROJECTION )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( MOMENTUM_PROJECTION_RHS )
KRATOS_DEFINE_VARIABLE( double, MASS_PROJECTION )
KRATOS_DEFINE_VARIABLE( double, MASS_PROJECTION_RHS )

}

#endif	/* KRATOS_STABILIZED_CFD_APPLICATION_VARIABLES_H_INCLUDED */
