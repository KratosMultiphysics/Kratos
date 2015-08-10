//
//   Project Name:        KratosAdjointFluidApplication $
//   Last modified by:    $Author: michael.andre@tum.de $
//   Date:                $Date:          February 2015 $
//   Revision:            $Revision:                0.0 $
//
//

#if !defined(KRATOS_ADJOINT_FLUID_APPLICATION_VARIABLES_H_INCLUDED)
#define  KRATOS_ADJOINT_FLUID_APPLICATION_VARIABLES_H_INCLUDED

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"

namespace Kratos
{
// Define variables
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( LAMBDA_VELOCITY )
KRATOS_DEFINE_VARIABLE(double, LAMBDA_PRESSURE )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( SHAPE_SENSITIVITY )
KRATOS_DEFINE_VARIABLE(double, NORMAL_SENSITIVITY )
}

#endif	/* KRATOS_ADJOINT_FLUID_APPLICATION_VARIABLES_H_INCLUDED */
