//
//   Project Name:        KratosParticleMechanicsApplication $
//   Created by:          $Author:                 IIaconeta $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:                   June 2017 $
//   Revision:            $Revision:                     0.0 $
//
//

// System includes

// External includes


// Project includes
#include "custom_constitutive/flow_rules/MPM_flow_rule.hpp"

namespace Kratos
{
 
    KRATOS_CREATE_LOCAL_FLAG( MPMFlowRule, IMPLEX_ACTIVE,                0 );
    KRATOS_CREATE_LOCAL_FLAG( MPMFlowRule, PLASTIC_REGION,               1 );
    KRATOS_CREATE_LOCAL_FLAG( MPMFlowRule, PLASTIC_RATE_REGION,          2 );
    KRATOS_CREATE_LOCAL_FLAG( MPMFlowRule, RETURN_MAPPING_COMPUTED,      3 );


}  // namespace Kratos.




