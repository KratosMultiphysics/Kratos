//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mohammad R. Hashemi
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
    // External interfacial force, e.g. for including the electromagentic coupling
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( DROPLET_DYNAMICS_APPLICATION, EXT_INT_FORCE)
}

#endif	/* KRATOS_DROPLET_DYNAMICS_APPLICATION_VARIABLES_H_INCLUDED */
