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
#include "custom_utilities/multipoint_constraint_data.hpp"

namespace Kratos
{
typedef MpcData::Pointer MpcDataPointerType;
typedef std::vector<MpcDataPointerType> *MpcDataPointerVectorType;

KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(CHIM_NEUMANN_COND);
// For MPC implementations
KRATOS_DEFINE_VARIABLE(MpcDataPointerVectorType, MPC_DATA_CONTAINER)
//KRATOS_DEFINE_VARIABLE(MpcDataPointerType, IS_WEAK);

// For segregated Chimera
KRATOS_DEFINE_VARIABLE(double, BOUNDARY_NODE);
KRATOS_DEFINE_VARIABLE(double, FLUX);
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(TRACTION);

}

#endif /* KRATOS_CHIMERA_APPLICATION_VARIABLES_H_INCLUDED */
