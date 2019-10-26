//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Tosi
//

#if !defined(KRATOS_EXAQUTE_SANDBOX_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_EXAQUTE_SANDBOX_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/kratos_application.h"

namespace Kratos
{
KRATOS_DEFINE_VARIABLE( double, DIVERGENCE_WEIGHTED )
KRATOS_DEFINE_VARIABLE( double, VELOCITY_H1_SEMINORM )
KRATOS_DEFINE_VARIABLE( double, PRESSURE_WEIGHTED )
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( EXAQUTE_SANDBOX_APPLICATION, VELOCITY_WEIGHTED )

}

#endif	/* KRATOS_EXAQUTE_SANDBOX_APPLICATION_VARIABLES_H_INCLUDED */
