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
template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
double CoulombFrictionalLaw<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::GetThresholdValue(
    const NodeType& rNode,
    const Condition& rCondition,
    const ProcessInfo& rCurrentProcessInfo
    )
{
//     // Manually compute the augmented normal pressure
//     const double common_epsilon = rCurrentProcessInfo[INITIAL_PENALTY];
//     const double scale_factor = rCurrentProcessInfo[SCALE_FACTOR];
//
//     const double epsilon = rNode.Has(INITIAL_PENALTY) ? rNode.GetValue(INITIAL_PENALTY) : common_epsilon;
//
//     const array_1d<double,3>& r_lagrange_multiplier = rNode.FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);
//     const array_1d<double,3>& r_nodal_normal = rNode.FastGetSolutionStepValue(NORMAL);
//     const double normal_r_lagrange_multiplier = inner_prod(r_nodal_normal, r_lagrange_multiplier);
//
//     const double augmented_normal_pressure = scale_factor * normal_r_lagrange_multiplier + epsilon * rNode.FastGetSolutionStepValue(WEIGHTED_GAP);

    // Getting the friction coefficient
    const double mu = this->GetFrictionCoefficient(rNode, rCondition, rCurrentProcessInfo);

    // Computing the tangent theshold pressure according to augmented normal pressure
    return - mu * rNode.GetValue(AUGMENTED_NORMAL_CONTACT_PRESSURE);
//     return - mu * augmented_normal_pressure;
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
double CoulombFrictionalLaw<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::GetDerivativeThresholdValue(
    const NodeType& rNode,
    const Condition& rCondition,
    const ProcessInfo& rCurrentProcessInfo,
    const DerivativeDataType& pDerivativeDatabase,
    const IndexType IndexDerivative
    )
{
    // Manually compute the augmented normal pressure
    const double common_epsilon = rCurrentProcessInfo[INITIAL_PENALTY];
    const double scale_factor = rCurrentProcessInfo[SCALE_FACTOR];

    const double epsilon = rNode.Has(INITIAL_PENALTY) ? rNode.GetValue(INITIAL_PENALTY) : common_epsilon;

    const array_1d<double,3>& r_lagrange_multiplier = rNode.FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);
    const array_1d<double,3>& r_nodal_normal = rNode.FastGetSolutionStepValue(NORMAL);
    const double normal_r_lagrange_multiplier = inner_prod(r_nodal_normal, r_lagrange_multiplier);

    const double augmented_normal_pressure = scale_factor * normal_r_lagrange_multiplier + epsilon * rNode.FastGetSolutionStepValue(WEIGHTED_GAP);

    const double mu = this->GetFrictionCoefficient(rNode, rCondition, rCurrentProcessInfo);

    const double derivative_threshold_value = - mu * (0.0); // TODO: Finnish this

    return derivative_threshold_value;
}

/***********************************************************************************/
/***********************************************************************************/

template class CoulombFrictionalLaw<2, 2, false, 2>;
template class CoulombFrictionalLaw<3, 3, false, 3>;
template class CoulombFrictionalLaw<3, 4, false, 4>;
template class CoulombFrictionalLaw<3, 3, false, 4>;
template class CoulombFrictionalLaw<3, 4, false, 3>;
template class CoulombFrictionalLaw<2, 2, true,  2>;
template class CoulombFrictionalLaw<3, 3, true,  3>;
template class CoulombFrictionalLaw<3, 4, true,  4>;
template class CoulombFrictionalLaw<3, 3, true,  4>;
template class CoulombFrictionalLaw<3, 4, true,  3>;

}  // namespace Kratos.

