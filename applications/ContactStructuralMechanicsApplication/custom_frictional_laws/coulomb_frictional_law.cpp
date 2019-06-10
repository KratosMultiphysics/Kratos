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
#include "custom_frictional_laws/coulomb_frictional_law.h"

namespace Kratos
{
double CoulombFrictionalLaw::GetThresholdValue(
    const NodeType& rNode,
    const Condition& rCondition,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const double mu = this->GetFrictionCoefficient(rNode, rCondition, rCurrentProcessInfo);
    return - mu * rNode.GetValue(AUGMENTED_NORMAL_CONTACT_PRESSURE);
}

}  // namespace Kratos.

