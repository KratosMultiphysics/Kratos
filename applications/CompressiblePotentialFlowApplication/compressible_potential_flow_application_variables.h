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
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, VELOCITY_INFINITY)
    KRATOS_DEFINE_VARIABLE( double, POSITIVE_POTENTIAL)
    KRATOS_DEFINE_VARIABLE( double, NEGATIVE_POTENTIAL)
    KRATOS_DEFINE_VARIABLE(Vector, WAKE_ELEMENTAL_DISTANCES )
    KRATOS_DEFINE_VARIABLE( double, WAKE_DISTANCE)
    KRATOS_DEFINE_VARIABLE(Vector, LEVEL_SET_ELEMENTAL_DISTANCES )
    KRATOS_DEFINE_VARIABLE( double, LEVEL_SET_DISTANCE)

}

#endif	/* KRATOS_COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION_VARIABLES_H_INCLUDED */
