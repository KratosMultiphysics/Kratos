//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Klaus B. Sautter
//

#if !defined(KRATOS_CABLE_NET_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_CABLE_NET_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "containers/variable.h"
#include "includes/ublas_interface.h"

namespace Kratos
{
    KRATOS_DEFINE_APPLICATION_VARIABLE( CABLE_NET_APPLICATION, Vector, SPRING_DEFORMATION_EMPIRICAL_POLYNOMIAL)
}

#endif	/* KRATOS_CABLE_NET_APPLICATION_VARIABLES_H_INCLUDED */
