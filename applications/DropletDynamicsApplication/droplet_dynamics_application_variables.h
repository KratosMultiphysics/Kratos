//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
//

#if !defined(KRATOS_DROPLET_DYNAMICS_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_DROPLET_DYNAMICS_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "containers/variable.h"

namespace Kratos
{
KRATOS_DEFINE_VARIABLE( double, DOF_1 )
KRATOS_DEFINE_VARIABLE( double, DOF_2 )
KRATOS_DEFINE_VARIABLE( double, ScalarVariable )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( VectorVariable )

}

#endif	/* KRATOS_DROPLET_DYNAMICS_APPLICATION_VARIABLES_H_INCLUDED */
