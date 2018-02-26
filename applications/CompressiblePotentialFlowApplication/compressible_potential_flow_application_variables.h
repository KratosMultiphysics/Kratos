//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(KRATOS_COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/kratos_application.h"

namespace Kratos
{

KRATOS_DEFINE_VARIABLE( bool, UPPER_SURFACE )
KRATOS_DEFINE_VARIABLE( bool, LOWER_SURFACE )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( WAKE_NORMAL )
KRATOS_DEFINE_VARIABLE( Matrix, PROJECTION_MATRIX )
KRATOS_DEFINE_VARIABLE( Matrix, UPPER_PROJECTION )
KRATOS_DEFINE_VARIABLE( Matrix, LOWER_PROJECTION )

}

#endif	/* KRATOS_COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION_VARIABLES_H_INCLUDED */
