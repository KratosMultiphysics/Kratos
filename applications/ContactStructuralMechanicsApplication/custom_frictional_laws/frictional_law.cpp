// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:         BSD License
//                   license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/variables.h"
#include "custom_frictional_laws/frictional_law.h"

namespace Kratos
{
double FrictionalLaw::GetFrictionCoefficient(
    const Node& rNode,
    const PairedCondition& rCondition,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const auto& r_properties = rCondition.GetProperties();
    if (r_properties.Has(FRICTION_COEFFICIENT)) {
        return r_properties.GetValue(FRICTION_COEFFICIENT);
    } else if (rNode.Has(FRICTION_COEFFICIENT)) {
        return rNode.GetValue(FRICTION_COEFFICIENT);
    } else if (rCurrentProcessInfo.Has(FRICTION_COEFFICIENT)) {
        return rCurrentProcessInfo.GetValue(FRICTION_COEFFICIENT);
    }

    return 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

double FrictionalLaw::GetThresholdValue(
    const Node& rNode,
    const PairedCondition& rCondition,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_ERROR << "You are calling to the base class method GetThresholdValue, check your frictional law declaration" << std::endl;

    return 0.0;
}

}  // namespace Kratos.

