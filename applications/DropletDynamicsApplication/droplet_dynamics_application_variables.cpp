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

#include "droplet_dynamics_application_variables.h"

namespace Kratos
{
    // External interfacial force, e.g. for including the electromagentic coupling
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(EXT_INT_FORCE)

    KRATOS_CREATE_VARIABLE( double, EPOTENTIAL )
    KRATOS_CREATE_VARIABLE( double, PERMITTIVITYPOS )
    KRATOS_CREATE_VARIABLE( double, PERMITTIVITYNEG )
    KRATOS_CREATE_VARIABLE( double, CONDUCTIVITYPOS )
    KRATOS_CREATE_VARIABLE( double, CONDUCTIVITYNEG )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( EFIELD )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( SCHARGE )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( EFORCE )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( EFIELDPOS )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( EFIELDNEG )
    // variables for enrich CASES
    KRATOS_CREATE_VARIABLE( double, INV_K_ENRICH )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( BIJ_ENRICH_ROW )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( BJI_ENRICH_ROW )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( POS_GRAD_ENRICH )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( NEG_GRAD_ENRICH )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( RHS_ENRICH )
}
