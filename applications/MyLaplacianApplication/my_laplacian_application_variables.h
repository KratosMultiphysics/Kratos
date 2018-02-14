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

#if !defined(KRATOS_MY_LAPLACIAN_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_MY_LAPLACIAN_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/kratos_application.h"

namespace Kratos
{
KRATOS_DEFINE_VARIABLE( double, MY_SCALAR )
KRATOS_DEFINE_VARIABLE( bool, MY_BOOL )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( MY_VECTOR )

}

#endif	/* KRATOS_MY_LAPLACIAN_APPLICATION_VARIABLES_H_INCLUDED */
