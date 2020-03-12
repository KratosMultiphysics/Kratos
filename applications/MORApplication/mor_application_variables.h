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

#if !defined(KRATOS_MOR_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_MOR_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/kratos_application.h"

namespace Kratos
{
KRATOS_DEFINE_VARIABLE( double, FREQUENCY )
KRATOS_DEFINE_VARIABLE( double, ACOUSTIC_PRESSURE )
KRATOS_DEFINE_VARIABLE( double, PRESSURE_GRADIENT )
KRATOS_DEFINE_VARIABLE( double, ACOUSTIC_VELOCITY )
KRATOS_DEFINE_VARIABLE( double, ACOUSTIC_PRESSURE_RESIDUAL)
// Nodal load variables
KRATOS_DEFINE_APPLICATION_VARIABLE(MOR_APPLICATION, double, ACOUSTIC_DISPLACEMENT )

}

#endif	/* KRATOS_MOR_APPLICATION_VARIABLES_H_INCLUDED */
