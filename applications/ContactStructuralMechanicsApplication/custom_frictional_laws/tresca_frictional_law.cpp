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
#include "contact_structural_mechanics_application_variables.h"
#include "custom_frictional_laws/tresca_frictional_law.h"

namespace Kratos
{
double TrescaFrictionalLaw::GetThresholdValue(
    const NodeType& rNode,
    const Condition& rCondition,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const auto& r_properties = rCondition.GetProperties();
    if (r_properties.Has(TRESCA_FRICTION_THRESHOLD)) {
        return r_properties.GetValue(TRESCA_FRICTION_THRESHOLD);
    } else if (rNode.Has(TRESCA_FRICTION_THRESHOLD)) {
        return rNode.GetValue(TRESCA_FRICTION_THRESHOLD);
    } else if (rCurrentProcessInfo.Has(TRESCA_FRICTION_THRESHOLD)) {
        return rCurrentProcessInfo.GetValue(TRESCA_FRICTION_THRESHOLD);
    }

    return 0.0;
}

}  // namespace Kratos.

