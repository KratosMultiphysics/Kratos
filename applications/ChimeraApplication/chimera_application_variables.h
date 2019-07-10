//
// ==============================================================================
//  ChimeraApplication
//
//  License:         BSD License
//                   license: ChimeraApplication/license.txt
//
//  Main authors:    Aditya Ghantasala, https://github.com/adityaghantasala
//                   Navaneeth K Narayanan
//
// ==============================================================================
//

#if !defined(KRATOS_CHIMERA_APPLICATION_VARIABLES_H_INCLUDED)
#define KRATOS_CHIMERA_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/kratos_application.h"
#include "chimera_application.h"

namespace Kratos
{

KRATOS_DEFINE_FLAG( FS_CHIMERA_VEL_CONSTRAINT);
KRATOS_DEFINE_FLAG( FS_CHIMERA_PRE_CONSTRAINT);
KRATOS_DEFINE_FLAG( CHIMERA_INTERNAL_BOUNDARY);
}

#endif /* KRATOS_CHIMERA_APPLICATION_VARIABLES_H_INCLUDED */
