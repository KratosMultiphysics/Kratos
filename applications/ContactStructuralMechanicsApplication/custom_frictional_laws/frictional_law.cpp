// KRATOS  ___|  |       |       |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//           | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License: BSD License
//   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/variables.h"
#include "custom_frictional_laws/frictional_law.h"

namespace Kratos
{
double FrictionalLaw::GetFrictionCoefficient(
    const NodeType& rNode,
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
    const NodeType& rNode,
    const PairedCondition& rCondition,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_ERROR << "You are calling to the base class method GetThresholdValue, check your frictional law declaration" << std::endl;

    return 0.0;
}

}  // namespace Kratos.

