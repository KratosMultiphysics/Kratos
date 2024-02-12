//
//   Project Name:        KratosPoromechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes


// Project includes
#include "custom_constitutive/custom_flow_rules/flow_rule.hpp"

namespace Kratos
{

    KRATOS_CREATE_LOCAL_FLAG( FlowRule, IMPLEX_ACTIVE,                0 );
    KRATOS_CREATE_LOCAL_FLAG( FlowRule, PLASTIC_REGION,               1 );
    KRATOS_CREATE_LOCAL_FLAG( FlowRule, PLASTIC_RATE_REGION,          2 );
    KRATOS_CREATE_LOCAL_FLAG( FlowRule, RETURN_MAPPING_COMPUTED,      3 );


}  // namespace Kratos.




