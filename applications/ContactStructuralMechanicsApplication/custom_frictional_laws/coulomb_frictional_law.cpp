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
#include "contact_structural_mechanics_application_variables.h"
#include "custom_frictional_laws/coulomb_frictional_law.h"

namespace Kratos
{
template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
double CoulombFrictionalLaw<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::GetThresholdValue(
    const Node& rNode,
    const PairedCondition& rCondition,
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
    const Node& rNode,
    const PairedCondition& rCondition,
    const ProcessInfo& rCurrentProcessInfo,
    const DerivativeDataType& rDerivativeData,
    const MortarConditionMatrices& rMortarConditionMatrices,
    const IndexType IndexDerivative,
    const IndexType IndexNode
    )
{
    const double mu = this->GetFrictionCoefficient(rNode, rCondition, rCurrentProcessInfo);

    double derivative_threshold_value = 0.0;

    const IndexType aux_index = IndexDerivative - TDim * (TNumNodes + TNumNodesMaster);
    if (aux_index > 0) {
        if (aux_index/TDim == IndexNode) {
            const double scale_factor = rCurrentProcessInfo[SCALE_FACTOR];
            const array_1d<double,3>& r_nodal_normal = rNode.FastGetSolutionStepValue(NORMAL);
            derivative_threshold_value = - mu * (scale_factor * r_nodal_normal[aux_index%TDim]);
        }
    } else {
        double delta_weighted_gap = 0.0;

        // Getting auxiliary indixes
        const IndexType derivative_node_index = IndexDerivative/TDim;
        const IndexType derivative_dimension_index = IndexDerivative%TDim;

        // Auxiliary values
        const double common_epsilon = rCurrentProcessInfo[INITIAL_PENALTY];
        const double epsilon = rNode.Has(INITIAL_PENALTY) ? rNode.GetValue(INITIAL_PENALTY) : common_epsilon;

        // Mortar condition matrices - DOperator and MOperator
        const BoundedMatrix<double, TNumNodes, TNumNodes>& DOperator = rMortarConditionMatrices.DOperator;
        const BoundedMatrix<double, TNumNodes, TNumNodesMaster>& MOperator = rMortarConditionMatrices.MOperator;
        const BoundedMatrix<double, TNumNodes, TNumNodes>& DeltaMOperator = rMortarConditionMatrices.DeltaMOperator[IndexDerivative];
        const BoundedMatrix<double, TNumNodes, TNumNodesMaster>& DeltaDOperator = rMortarConditionMatrices.DeltaDOperator[IndexDerivative];

        // Some operations
        array_1d<double, TDim> aux_normal, aux_array;
        const array_1d<double,3>& r_nodal_normal = rNode.FastGetSolutionStepValue(NORMAL);

        for (IndexType i_dim = 0; i_dim < TDim; ++i_dim) {
            aux_normal[i_dim] = r_nodal_normal[i_dim];
        }

        // Delta coordinates
        BoundedMatrix<double, TNumNodes, TDim> Deltax1 = ZeroMatrix(TNumNodes, TDim);
        BoundedMatrix<double, TNumNodesMaster, TDim> Deltax2 = ZeroMatrix(TNumNodesMaster, TDim);

        // Update matrix delta coordinates
        if (derivative_node_index < TNumNodes) {
            Deltax1(derivative_node_index, derivative_dimension_index) = 1.0;
        } else {
            Deltax2(derivative_node_index - TNumNodes, derivative_dimension_index) = 1.0;
        }

        // Adding delta contribution
        const BoundedMatrix<double, TNumNodes, TDim> D_Deltax1_M_Deltax2 = prod(MOperator, Deltax1) - prod(DOperator, Deltax2);
        noalias(aux_array) = row(D_Deltax1_M_Deltax2, IndexNode);
        delta_weighted_gap += inner_prod(aux_array, - aux_normal);

        // Current coordinates
        const BoundedMatrix<double, TNumNodes, TDim> x1 = rDerivativeData.u1 + rDerivativeData.X1;
        const BoundedMatrix<double, TNumNodesMaster, TDim> x2 = rDerivativeData.u2 + rDerivativeData.X2;

        // Adding delta contribution
        const BoundedMatrix<double, TNumNodes, TDim> DeltaD_x1_DeltaM_x2 = prod(DeltaMOperator, x1) - prod(DeltaDOperator, x2);
        noalias(aux_array) = row(DeltaD_x1_DeltaM_x2, IndexNode);
        delta_weighted_gap += inner_prod(aux_array, - aux_normal);

        // Delta normal contribution
        if (TNormalVariation && (derivative_node_index < TNumNodes)) {
            const auto& DeltaNormalSlave = rDerivativeData.DeltaNormalSlave[IndexDerivative];
            noalias(aux_normal) = row(DeltaNormalSlave, IndexNode);
            const BoundedMatrix<double, TNumNodes, TDim> D_x1_M_x2 = prod(MOperator, x1) - prod(DOperator, x2);
            noalias(aux_array) = row(D_x1_M_x2, IndexNode);
            delta_weighted_gap += inner_prod(aux_array, - aux_normal);
        }

        derivative_threshold_value = - mu * (epsilon * delta_weighted_gap);
    }

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

