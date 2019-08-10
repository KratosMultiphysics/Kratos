//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:   Ilaria Iaconeta
//


// System includes

// External includes

// Project includes
#include "custom_constitutive/flow_rules/particle_flow_rule.hpp"

namespace Kratos
{

KRATOS_CREATE_LOCAL_FLAG( ParticleFlowRule, IMPLEX_ACTIVE,                0 );
KRATOS_CREATE_LOCAL_FLAG( ParticleFlowRule, PLASTIC_REGION,               1 );
KRATOS_CREATE_LOCAL_FLAG( ParticleFlowRule, PLASTIC_RATE_REGION,          2 );
KRATOS_CREATE_LOCAL_FLAG( ParticleFlowRule, RETURN_MAPPING_COMPUTED,      3 );


}  // namespace Kratos.




