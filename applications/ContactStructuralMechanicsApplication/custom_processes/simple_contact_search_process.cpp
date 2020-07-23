// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
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
    typename NodesArrayType::iterator ItNode,
    const double CommonEpsilon,
    const double ScaleFactor
    )
{
    // First we activate
    BaseType::SetActiveNode(ItNode, CommonEpsilon);

    // Normal gap
    const double normal_gap = ItNode->Has(NORMAL_GAP) ? ItNode->GetValue(NORMAL_GAP) : 0.0;

    // In case of penetration
    if (normal_gap < 0.0) {
        // Auxiliar values
        const double epsilon = (ItNode->Has(INITIAL_PENALTY) ? ItNode->GetValue(INITIAL_PENALTY) : CommonEpsilon)/ScaleFactor;
        const double nodal_area = ItNode->Has(NODAL_AREA) ? ItNode->GetValue(NODAL_AREA) : 1.0;

        // Setting approximation
        switch(BaseType::mTypeSolution) {
            case BaseType::TypeSolution::VectorLagrangeMultiplier :
                noalias(ItNode->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER)) = epsilon * nodal_area * normal_gap * ItNode->FastGetSolutionStepValue(NORMAL);
                break;
            case BaseType::TypeSolution::ScalarLagrangeMultiplier :
                ItNode->FastGetSolutionStepValue(SCALAR_LAGRANGE_MULTIPLIER) = epsilon * nodal_area * normal_gap;
                break;
            case BaseType::TypeSolution::NormalContactStress :
                ItNode->FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE) = epsilon * nodal_area * normal_gap;
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
