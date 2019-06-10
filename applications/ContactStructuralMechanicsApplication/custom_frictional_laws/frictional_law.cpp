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
template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
double FrictionalLaw<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::GetFrictionCoefficient(
    const NodeType& rNode,
    const Condition& rCondition,
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

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
double FrictionalLaw<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::GetThresholdValue(
    const NodeType& rNode,
    const Condition& rCondition,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_ERROR << "You are calling to the base class method GetThresholdValue, check your frictional law declaration" << std::endl;

    return 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
double FrictionalLaw<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::GetDerivativeThresholdValue(
    const NodeType& rNode,
    const Condition& rCondition,
    const ProcessInfo& rCurrentProcessInfo,
    const DerivativeDataType& pDerivativeDatabase,
    const IndexType IndexDerivative
    )
{
    KRATOS_ERROR << "You are calling to the base class method GetDerivativeThresholdValue, check your frictional law declaration" << std::endl;

    return 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

template class FrictionalLaw<2, 2, false, 2>;
template class FrictionalLaw<3, 3, false, 3>;
template class FrictionalLaw<3, 4, false, 4>;
template class FrictionalLaw<3, 3, false, 4>;
template class FrictionalLaw<3, 4, false, 3>;
template class FrictionalLaw<2, 2, true,  2>;
template class FrictionalLaw<3, 3, true,  3>;
template class FrictionalLaw<3, 4, true,  4>;
template class FrictionalLaw<3, 3, true,  4>;
template class FrictionalLaw<3, 4, true,  3>;

}  // namespace Kratos.

