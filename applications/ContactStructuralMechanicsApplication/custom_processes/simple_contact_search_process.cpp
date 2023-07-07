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
#include "utilities/variable_utils.h"
#include "custom_processes/simple_contact_search_process.h"
#include "contact_structural_mechanics_application_variables.h"

namespace Kratos
{
template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
SimpleContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>::SimpleContactSearchProcess(
    ModelPart & rMainModelPart,
    Parameters ThisParameters,
    Properties::Pointer pPairedProperties
    ) : BaseType(rMainModelPart, ThisParameters, pPairedProperties)
{
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void SimpleContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>::SetActiveNode(
    Node& rNode,
    const double CommonEpsilon,
    const double ScaleFactor
    )
{
    // First we activate
    BaseType::SetActiveNode(rNode, CommonEpsilon);

    // Normal gap
    const double normal_gap = rNode.Has(NORMAL_GAP) ? rNode.GetValue(NORMAL_GAP) : 0.0;

    // In case of penetration
    if (normal_gap < 0.0) {
        // Auxiliary values
        const double epsilon = (rNode.Has(INITIAL_PENALTY) ? rNode.GetValue(INITIAL_PENALTY) : CommonEpsilon)/ScaleFactor;
        const double nodal_area = rNode.Has(NODAL_AREA) ? rNode.GetValue(NODAL_AREA) : 1.0;

        // Setting approximation
        switch(BaseType::mTypeSolution) {
            case BaseType::TypeSolution::VectorLagrangeMultiplier :
                noalias(rNode.FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER)) = epsilon * nodal_area * normal_gap * rNode.FastGetSolutionStepValue(NORMAL);
                break;
            case BaseType::TypeSolution::ScalarLagrangeMultiplier :
                rNode.FastGetSolutionStepValue(SCALAR_LAGRANGE_MULTIPLIER) = epsilon * nodal_area * normal_gap;
                break;
            case BaseType::TypeSolution::NormalContactStress :
                rNode.FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE) = epsilon * nodal_area * normal_gap;
                break;
            case BaseType::TypeSolution::FrictionlessPenaltyMethod :
                break;
            case BaseType::TypeSolution::FrictionalPenaltyMethod :
                break;
            case BaseType::TypeSolution::OtherFrictionless :
                break;
            case BaseType::TypeSolution::OtherFrictional :
                break;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template class SimpleContactSearchProcess<2, 2>;
template class SimpleContactSearchProcess<3, 3>;
template class SimpleContactSearchProcess<3, 4>;
template class SimpleContactSearchProcess<3, 3, 4>;
template class SimpleContactSearchProcess<3, 4, 3>;

}  // namespace Kratos.
